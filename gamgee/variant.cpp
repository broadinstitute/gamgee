#include "variant.h"
#include "individual_field.h"
#include "individual_field_value.h"
#include "shared_field.h"
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

IndividualField<IndividualFieldValue<int32_t>> Variant::genotype_quals() const {
  return integer_individual_field("GQ");
}

IndividualField<IndividualFieldValue<int32_t>> Variant::phred_likelihoods() const {
  return integer_individual_field("PL");
}

IndividualField<IndividualFieldValue<int32_t>> Variant::integer_individual_field(const std::string& tag) const {
  const auto id = bcf_hdr_id2int(m_header.get(), BCF_DT_ID, tag.c_str());
  if (!bcf_hdr_idinfo_exists(m_header.get(), BCF_HL_FMT, id))
    return IndividualField<IndividualFieldValue<int32_t>>{}; 
  if (bcf_hdr_id2type(m_header.get(),BCF_HL_FMT,id)!=BCF_HT_INT)
    throw runtime_error("individual field requested is not an integer");
  return individual_field_as_integer(tag); // @todo: move this to an indexed based lookup API
}

IndividualField<IndividualFieldValue<float>> Variant::float_individual_field(const std::string& tag) const {
  const auto id = bcf_hdr_id2int(m_header.get(), BCF_DT_ID, tag.c_str());
  if (!bcf_hdr_idinfo_exists(m_header.get(), BCF_HL_FMT, id))
    return IndividualField<IndividualFieldValue<float>>{};
  if (bcf_hdr_id2type(m_header.get(),BCF_HL_FMT,id)!=BCF_HT_REAL)
    throw runtime_error("individual field requested is not a float");
  return individual_field_as_float(tag); // @todo: move this to an indexed based lookup API
}

IndividualField<IndividualFieldValue<string>> Variant::string_individual_field(const std::string& tag) const {
  const auto id = bcf_hdr_id2int(m_header.get(), BCF_DT_ID, tag.c_str());
  if (!bcf_hdr_idinfo_exists(m_header.get(), BCF_HL_FMT, id))
    return IndividualField<IndividualFieldValue<string>>{};
  if (bcf_hdr_id2type(m_header.get(),BCF_HL_FMT,id)!=BCF_HT_STR)
    throw runtime_error("individual field requested is not a string");
  return individual_field_as_string(tag); // @todo: move this to an indexed based lookup API
}

IndividualField<IndividualFieldValue<int32_t>> Variant::individual_field_as_integer(const std::string& tag) const {
  const auto fmt = find_individual_field_by_tag(tag);
  if (fmt == nullptr) ///< if the variant is missing or the PL tag is missing, return an empty IndividualField
    return IndividualField<IndividualFieldValue<int32_t>>{};
  return IndividualField<IndividualFieldValue<int32_t>>{m_body, fmt};
}

IndividualField<IndividualFieldValue<float>> Variant::individual_field_as_float(const std::string& tag) const {
  const auto fmt = find_individual_field_by_tag(tag);
  if (fmt == nullptr) ///< if the variant is missing or the PL tag is missing, return an empty IndividualField
    return IndividualField<IndividualFieldValue<float>>{};
  return IndividualField<IndividualFieldValue<float>>{m_body, fmt};
}

IndividualField<IndividualFieldValue<std::string>> Variant::individual_field_as_string(const std::string& tag) const {
  const auto fmt = find_individual_field_by_tag(tag);
  if (fmt == nullptr) ///< if the variant is missing or the PL tag is missing, return an empty IndividualField
    return IndividualField<IndividualFieldValue<string>>{};
  return IndividualField<IndividualFieldValue<string>>{m_body, fmt};
}

inline bcf_fmt_t* Variant::find_individual_field_by_tag(const string& tag) const {
  return bcf_get_fmt(m_header.get(), m_body.get(), tag.c_str());     
}

inline bcf_info_t* Variant::find_shared_field_by_tag(const string& tag) const {
  return bcf_get_info(m_header.get(), m_body.get(), tag.c_str());     
}

bool Variant::boolean_shared_field(const std::string& tag) const {
  const auto info = find_shared_field_by_tag(tag);
  return info != nullptr;
}

SharedField<int32_t> Variant::integer_shared_field(const std::string& tag) const {
  const auto id = bcf_hdr_id2int(m_header.get(), BCF_DT_ID, tag.c_str());
  if (!bcf_hdr_idinfo_exists(m_header.get(), BCF_HL_INFO, id))
    return SharedField<int32_t>{};
  if (bcf_hdr_id2type(m_header.get(), BCF_HL_INFO,id) != BCF_HT_INT)
    throw runtime_error("shared field requested is not a int");
  return shared_field_as_integer(tag); // @todo: move this to an indexed based lookup API
}

SharedField<float> Variant::float_shared_field(const std::string& tag) const {
  const auto id = bcf_hdr_id2int(m_header.get(), BCF_DT_ID, tag.c_str());
  if (!bcf_hdr_idinfo_exists(m_header.get(), BCF_HL_INFO, id))
    return SharedField<float>{};
  if (bcf_hdr_id2type(m_header.get(), BCF_HL_INFO,id) != BCF_HT_REAL)
    throw runtime_error("shared field requested is not a float");
  return shared_field_as_float(tag); // @todo: move this to an indexed based lookup API
}

SharedField<string> Variant::string_shared_field(const std::string& tag) const {
  const auto id = bcf_hdr_id2int(m_header.get(), BCF_DT_ID, tag.c_str());
  if (!bcf_hdr_idinfo_exists(m_header.get(), BCF_HL_INFO, id)) 
    return SharedField<string>{};
  if (bcf_hdr_id2type(m_header.get(), BCF_HL_INFO,id) != BCF_HT_STR)
    throw runtime_error("shared field requested is not a string");
  return shared_field_as_string(tag); // @todo: move this to an indexed based lookup API
}

SharedField<int32_t> Variant::shared_field_as_integer(const std::string& tag) const {
  const auto info = find_shared_field_by_tag(tag);
  if (info == nullptr) 
    return SharedField<int32_t>{};
  return SharedField<int32_t>{m_body, info};
}

SharedField<float> Variant::shared_field_as_float(const std::string& tag) const {
  const auto info = find_shared_field_by_tag(tag);
  if (info == nullptr)
    return SharedField<float>{};
  return SharedField<float>{m_body, info};
}

SharedField<string> Variant::shared_field_as_string(const std::string& tag) const {
  const auto info = find_shared_field_by_tag(tag);
  if (info == nullptr)
    return SharedField<string>{};
  return SharedField<string>{m_body, info};
}

IndividualField<Genotype> Variant::genotypes() const {
  // bcf_get_fmt() will unpack the record if necessary
  const auto fmt = bcf_get_fmt(m_header.get(), m_body.get(), "GT");
  if (fmt == nullptr) ///< if the variant is missing or the GT tag is missing, return an empty IndividualField
    return IndividualField<Genotype>{};
  return IndividualField<Genotype>{m_body, fmt};
}

}
