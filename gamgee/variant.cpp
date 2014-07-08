#include "variant.h"
#include "variant_field.h"
#include "variant_field_value.h"
#include "genotype_field_value.h"
#include "utils/hts_memory.h"
#include "utils/utils.h"

#include "htslib/vcf.h"

using namespace std;
namespace gamgee {

/**
 * @brief creates a variant record that points to htslib memory already allocated
 * @note the resulting Variant shares ownership of the pre-allocated memory via shared_ptr reference counting
 */
Variant::Variant(const std::shared_ptr<bcf_hdr_t>& header, const std::shared_ptr<bcf1_t>& body) noexcept :
  m_header {header},
  m_body {body}
{}

/**
 * @brief creates a deep copy of a variant record and header
 */
Variant::Variant(const Variant& other) :
  m_header {utils::make_shared_variant_header(utils::variant_header_deep_copy(other.m_header.get()))},
  m_body {utils::make_shared_variant(utils::variant_deep_copy(other.m_body.get()))}
{}

/**
 * @brief moves a variant record and header, transferring ownership of the underlying htslib memory
 */
Variant::Variant(Variant&& other) noexcept :
  m_header {move(other.m_header)},
  m_body {move(other.m_body)}
{}

/**
 * @brief creates a deep copy of a variant record
 * @param other the Variant to be copied
 * @note the parameter is passed by copy so we can maintain the strong exception safety guarantee
 */
Variant& Variant::operator=(const Variant& other) {
  if ( &other == this )  
    return *this;
  m_body = utils::make_shared_variant(utils::variant_deep_copy(other.m_body.get()));                   ///< shared_ptr assignment will take care of deallocating old record if necessary
  m_header = utils::make_shared_variant_header(utils::variant_header_deep_copy(other.m_header.get())); ///< shared_ptr assignment will take care of deallocating old record if necessary
  return *this;
}

/**
 * @brief moves a variant record, transferring ownership of the underlying htslib memory
 */
Variant& Variant::operator=(Variant&& other) noexcept {
  if ( &other == this )  
    return *this;
  m_body = move(other.m_body);
  m_header = move(other.m_header);
  return *this;
}

std::string Variant::ref() const {
  bcf_unpack(m_body.get(), BCF_UN_STR);
  return n_alleles() > 0 ?  string{m_body.get()->d.allele[0]} : string{}; 
}

std::vector<std::string> Variant::alt() const {
  bcf_unpack(m_body.get(), BCF_UN_STR);
  const auto n_all = n_alleles();
  return n_all > 1 ? utils::hts_string_array_to_vector(m_body.get()->d.allele+1, n_all-1) : vector<string>{}; // skip the first allele because it's the ref
}

VariantFilters Variant::filters() const {
  return VariantFilters{m_header, m_body};
}

bool Variant::has_filter(const std::string& filter) const {
  return bcf_has_filter(m_header.get(), m_body.get(), const_cast<char*>(filter.c_str())) > 0; // have to cast away the constness here for the C api to work. But the promise still remains as the C function is not modifying the string.
}

bool Variant::is_this_genotype(const DiploidPLGenotype& genotype, const uint32_t sample_index) const {
  const auto pl = phred_likelihoods();
  return !pl.empty() && pl[sample_index][static_cast<int32_t>(genotype)] == 0;
}

bool Variant::is_hom_ref(const uint32_t sample_index) const {
  return is_this_genotype(DiploidPLGenotype::HOM_REF, sample_index);
}

bool Variant::is_het(const uint32_t sample_index) const {
  return is_this_genotype(DiploidPLGenotype::HET, sample_index);
}

bool Variant::is_hom_var(const uint32_t sample_index) const {
  return is_this_genotype(DiploidPLGenotype::HOM_VAR, sample_index);
}

VariantField<VariantFieldValue<int32_t>> Variant::genotype_quals() const {
  return generic_integer_format_field("GQ");
}

VariantField<VariantFieldValue<int32_t>> Variant::phred_likelihoods() const {
  return generic_integer_format_field("PL");
}

VariantField<VariantFieldValue<int32_t>> Variant::generic_integer_format_field(const std::string& tag) const {
  const auto fmt = find_format_field_by_tag(tag);
  if (fmt == nullptr) ///< if the variant is missing or the PL tag is missing, return an empty VariantField
    return VariantField<VariantFieldValue<int32_t>>{};
  return VariantField<VariantFieldValue<int32_t>>{m_body, fmt};
}

VariantField<VariantFieldValue<float>> Variant::generic_float_format_field(const std::string& tag) const {
  const auto fmt = find_format_field_by_tag(tag);
  if (fmt == nullptr) ///< if the variant is missing or the PL tag is missing, return an empty VariantField
    return VariantField<VariantFieldValue<float>>{};
  return VariantField<VariantFieldValue<float>>{m_body, fmt};
}

VariantField<VariantFieldValue<std::string>> Variant::generic_string_format_field(const std::string& tag) const {
  const auto fmt = find_format_field_by_tag(tag);
  if (fmt == nullptr) ///< if the variant is missing or the PL tag is missing, return an empty VariantField
    return VariantField<VariantFieldValue<string>>{};
  return VariantField<VariantFieldValue<string>>{m_body, fmt};
}

inline bcf_fmt_t* Variant::find_format_field_by_tag(const string& tag) const {
  return bcf_get_fmt(m_header.get(), m_body.get(), tag.c_str());     
}

std::vector<int32_t> Variant::generic_integer_info_field(const std::string& tag) const {
  return generic_info_field<int32_t>(tag, BCF_HT_INT);
}

std::vector<float> Variant::generic_float_info_field(const std::string& tag) const {
  return generic_info_field<float>(tag, BCF_HT_REAL);
}

template <typename TYPE> inline std::vector<TYPE> Variant::generic_info_field(const std::string& tag, const int type) const {
  // Using malloc instead of new since bcf_get_info_values does realloc: http://www.stroustrup.com/bs_faq2.html#realloc
  auto mem = (TYPE *)malloc(sizeof(TYPE)); // bcf_get_info_values writes without realloc when count returns with 1;
  auto count = 1;
  const auto info_result = bcf_get_info_values(m_header.get(), m_body.get(), tag.c_str(), (void**)&mem, &count, type);
  if (info_result < 0) {
    return std::vector<TYPE>{}; // return empty for all errors, even asking for the wrong type
  }
  const auto results = std::vector<TYPE>(mem, mem + count); // int32_t and floats returned as arrays
  free(mem);
  return results;
}

std::vector<std::string> Variant::generic_string_info_field(const std::string& tag) const {
  // Using malloc instead of new since bcf_get_info_values does realloc: http://www.stroustrup.com/bs_faq2.html#realloc
  auto mem = (char*)malloc(0);
  auto count = 0;
  const auto info_result = bcf_get_info_string(m_header.get(), m_body.get(), tag.c_str(), &mem, &count);
  if (info_result < 0) {
    return std::vector<string>{}; // return empty for all errors, even asking for the wrong type
  }
  const auto results = std::vector<string>{std::string{mem}}; // strings returned as a single string
  free(mem);
  return results;
}

std::vector<bool> Variant::generic_boolean_info_field(const std::string& tag) const {
  const auto info_result = bcf_get_info_flag(m_header.get(), m_body.get(), tag.c_str(), NULL, NULL);
  if (info_result < 0) {
    return std::vector<bool>{}; // return empty for all errors, even asking for the wrong type
  }
  const auto results = std::vector<bool>{info_result == 1}; // flags are returned only in the return value
  return results;
}

VariantField<GenotypeFieldValue> Variant::phred_likelihoods_field() const {
  const auto fmt = bcf_get_fmt(m_header.get(), m_body.get(), "PL");
  if (fmt == nullptr) ///< if the variant is missing or the PL tag is missing, return an empty VariantField
    return VariantField<GenotypeFieldValue>{};
  return VariantField<GenotypeFieldValue>{m_body, fmt};
}

}
