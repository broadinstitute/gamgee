#include "variant.h"
#include "variant_field.h"
#include "variant_field_value.h"
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
 * @brief creates a deep copy of a variant record
 *
 * @note does not perform a deep copy of the variant header; to copy the header,
 *       first get it via the header() function and then copy it via the usual C++
 *       semantics
 */
Variant::Variant(const Variant& other) :
  m_header {other.m_header},
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
 * @note does not perform a deep copy of the variant header; to copy the header,
 *       first get it via the header() function and then copy it via the usual C++
 *       semantics
 */
Variant& Variant::operator=(const Variant& other) {
  if ( &other == this )  
    return *this;
  m_header = other.m_header;    ///< shared_ptr assignment will take care of deallocating old record if necessary
  m_body = utils::make_shared_variant(utils::variant_deep_copy(other.m_body.get()));  ///< shared_ptr assignment will take care of deallocating old record if necessary
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

std::string Variant::id () const {
  bcf_unpack(m_body.get(), BCF_UN_STR);
  return std::string{m_body->d.id};
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
  bcf_unpack(m_body.get(), BCF_UN_FLT);
  return VariantFilters{m_header, m_body};
}

bool Variant::has_filter(const std::string& filter) const {
  return bcf_has_filter(m_header.get(), m_body.get(), const_cast<char*>(filter.c_str())) > 0; // have to cast away the constness here for the C api to work. But the promise still remains as the C function is not modifying the string.
}

VariantField<VariantFieldValue<int32_t>> Variant::genotype_quals() const {
  return integer_individual_field("GQ");
}

VariantField<VariantFieldValue<int32_t>> Variant::phred_likelihoods() const {
  return integer_individual_field("PL");
}

VariantField<VariantFieldValue<int32_t>> Variant::integer_individual_field(const std::string& tag) const {
  const auto fmt = find_individual_field_by_tag(tag);
  if (fmt == nullptr) ///< if the variant is missing or the PL tag is missing, return an empty VariantField
    return VariantField<VariantFieldValue<int32_t>>{};
  return VariantField<VariantFieldValue<int32_t>>{m_body, fmt};
}

VariantField<VariantFieldValue<float>> Variant::float_individual_field(const std::string& tag) const {
  const auto fmt = find_individual_field_by_tag(tag);
  if (fmt == nullptr) ///< if the variant is missing or the PL tag is missing, return an empty VariantField
    return VariantField<VariantFieldValue<float>>{};
  return VariantField<VariantFieldValue<float>>{m_body, fmt};
}

VariantField<VariantFieldValue<std::string>> Variant::string_individual_field(const std::string& tag) const {
  const auto fmt = find_individual_field_by_tag(tag);
  if (fmt == nullptr) ///< if the variant is missing or the PL tag is missing, return an empty VariantField
    return VariantField<VariantFieldValue<string>>{};
  return VariantField<VariantFieldValue<string>>{m_body, fmt};
}

inline bcf_fmt_t* Variant::find_individual_field_by_tag(const string& tag) const {
  return bcf_get_fmt(m_header.get(), m_body.get(), tag.c_str());     
}

std::vector<int32_t> Variant::integer_shared_field(const std::string& tag) const {
  return shared_field<int32_t>(tag, BCF_HT_INT);
}

std::vector<float> Variant::float_shared_field(const std::string& tag) const {
  return shared_field<float>(tag, BCF_HT_REAL);
}

template <typename TYPE> 
inline std::vector<TYPE> Variant::shared_field(const std::string& tag, const int type) const {
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

std::vector<std::string> Variant::string_shared_field(const std::string& tag) const {
  // Using malloc instead of new since bcf_get_info_values does realloc: http://www.stroustrup.com/bs_faq2.html#realloc
  auto mem = (char*)malloc(0);
  auto count = 0;
  const auto info_result = bcf_get_info_string(m_header.get(), m_body.get(), tag.c_str(), &mem, &count);
  if (info_result < 0) {
    free(mem);
    return std::vector<string>{}; // return empty for all errors, even asking for the wrong type
  }
  const auto results = std::vector<string>{std::string{mem}}; // strings returned as a single string
  free(mem);
  return results;
}

std::vector<bool> Variant::boolean_shared_field(const std::string& tag) const {
  const auto info_result = bcf_get_info_flag(m_header.get(), m_body.get(), tag.c_str(), NULL, NULL);
  if (info_result < 0) {
    return std::vector<bool>{}; // return empty for all errors, even asking for the wrong type
  }
  const auto results = std::vector<bool>{info_result == 1}; // flags are returned only in the return value
  return results;
}

VariantField<Genotype> Variant::genotypes() const {
  // bcf_get_fmt() will unpack the record if necessary
  const auto fmt = bcf_get_fmt(m_header.get(), m_body.get(), "GT");
  if (fmt == nullptr) ///< if the variant is missing or the GT tag is missing, return an empty VariantField
    return VariantField<Genotype>{};
  return VariantField<Genotype>{m_body, fmt};
}

}
