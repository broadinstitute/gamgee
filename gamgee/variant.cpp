#include "variant.h"
#include "individual_field.h"
#include "individual_field_value.h"
#include "shared_field.h"
#include "utils/hts_memory.h"
#include "utils/utils.h"

#include "htslib/vcf.h"

using namespace std;
namespace gamgee {

/******************************************************************************
 * Accessory functions (private functions)                                    *
 ******************************************************************************/
inline bool Variant::check_field(const int32_t type_field, const int32_t type_value, const int32_t index) const {
  if ( ! m_header.has_field(index, type_field) )
    return false;
  if ( m_header.field_type(index, type_field) != type_value )
    throw std::runtime_error("field requested is not of the right type");
  return true;
}

template<class FIELD_TYPE, class INDEX_OR_TAG>
IndividualField<IndividualFieldValue<FIELD_TYPE>> Variant::individual_field_as(const INDEX_OR_TAG& p) const {
  const auto field_ptr = find_individual_field(p);
  if (field_ptr == nullptr)
    return IndividualField<IndividualFieldValue<FIELD_TYPE>>{};
  return IndividualField<IndividualFieldValue<FIELD_TYPE>>{m_body, field_ptr};
}

template<class FIELD_TYPE, class INDEX_OR_TAG>
SharedField<FIELD_TYPE> Variant::shared_field_as(const INDEX_OR_TAG& p) const {
  const auto field_ptr = find_shared_field(p);
  if (field_ptr == nullptr)
    return SharedField<FIELD_TYPE>{};
  return SharedField<FIELD_TYPE>{m_body, field_ptr};
}

AlleleType Variant::allele_type_from_difference(const int diff) const {
  if (diff == 0)
    return AlleleType::SNP;
  return (diff > 0) ? AlleleType::INSERTION : AlleleType::DELETION;
}



/******************************************************************************
 * Constructors and operator overloads                                        *
 ******************************************************************************/
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
  m_header {other.m_header.m_header},   // Avoid a deep copy here by constructing using other's internal shared header pointer
  m_body {utils::make_shared_variant(utils::variant_deep_copy(other.m_body.get()))}
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
  m_header = VariantHeader{other.m_header.m_header};    // Avoid a deep copy here by constructing using other's internal shared header pointer
  m_body = utils::make_shared_variant(utils::variant_deep_copy(other.m_body.get()));  ///< shared_ptr assignment will take care of deallocating old record if necessary
  return *this;
}


/******************************************************************************
 * General record API                                                         *
 ******************************************************************************/
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
  return VariantFilters{m_header.m_header, m_body};
}

bool Variant::has_filter(const std::string& filter) const {
  return bcf_has_filter(m_header.m_header.get(), m_body.get(), const_cast<char*>(filter.c_str())) > 0; // have to cast away the constness here for the C api to work. But the promise still remains as the C function is not modifying the string.
}

AlleleMask Variant::allele_mask() const {
  auto result = AlleleMask{};
  result.reserve(n_alleles());
  result.emplace_back(AlleleType::REFERENCE); // add the reference first
  const auto ref_allele = ref();
  for (const auto& alt_allele : alt()) {
    const auto diff = int32_t(alt_allele.size() - ref_allele.size());
    result.emplace_back(allele_type_from_difference(diff));
  }
  return result;
}

/******************************************************************************
 * Individual field API                                                       *
 ******************************************************************************/
IndividualField<IndividualFieldValue<int32_t>> Variant::integer_individual_field(const std::string& tag) const {
  return integer_individual_field(m_header.field_index(tag));
}

IndividualField<IndividualFieldValue<float>> Variant::float_individual_field(const std::string& tag) const {
  return float_individual_field(m_header.field_index(tag));
}

IndividualField<IndividualFieldValue<string>> Variant::string_individual_field(const std::string& tag) const {
  return string_individual_field(m_header.field_index(tag));
}

IndividualField<IndividualFieldValue<int32_t>> Variant::individual_field_as_integer(const std::string& tag) const {
  return individual_field_as<int32_t>(tag);
}

IndividualField<IndividualFieldValue<float>> Variant::individual_field_as_float(const std::string& tag) const {
  return individual_field_as<float>(tag);
}

IndividualField<IndividualFieldValue<std::string>> Variant::individual_field_as_string(const std::string& tag) const {
  return individual_field_as<string>(tag);
}

IndividualField<IndividualFieldValue<int32_t>> Variant::integer_individual_field(const int32_t index) const {
  if (check_field(BCF_HL_FMT, BCF_HT_INT, index))
    return individual_field_as<int32_t>(index); 
  return IndividualField<IndividualFieldValue<int32_t>>{};
}

IndividualField<IndividualFieldValue<float>> Variant::float_individual_field(const int32_t index) const {
  if (check_field(BCF_HL_FMT, BCF_HT_REAL, index))
    return individual_field_as<float>(index); 
  return IndividualField<IndividualFieldValue<float>>{};
}

IndividualField<IndividualFieldValue<string>> Variant::string_individual_field(const int32_t index) const {
  if (check_field(BCF_HL_FMT, BCF_HT_STR, index))
    return individual_field_as<string>(index); 
  return IndividualField<IndividualFieldValue<string>>{};
}

IndividualField<IndividualFieldValue<int32_t>> Variant::individual_field_as_integer(const int32_t index) const {
  return individual_field_as<int32_t>(index);
}

IndividualField<IndividualFieldValue<float>> Variant::individual_field_as_float(const int32_t index) const {
  return individual_field_as<float>(index);
}

IndividualField<IndividualFieldValue<std::string>> Variant::individual_field_as_string(const int32_t index) const {
  return individual_field_as<string>(index);
}


/******************************************************************************
 * Shared field API                                                           *
 ******************************************************************************/
bool Variant::boolean_shared_field(const std::string& tag) const {
  return find_shared_field(tag) != nullptr;
}

bool Variant::boolean_shared_field(const int32_t index) const {
  return find_shared_field(index) != nullptr;
}

SharedField<int32_t> Variant::integer_shared_field(const std::string& tag) const {
  return integer_shared_field(m_header.field_index(tag));
}

SharedField<float> Variant::float_shared_field(const std::string& tag) const {
  return float_shared_field(m_header.field_index(tag));
}

SharedField<string> Variant::string_shared_field(const std::string& tag) const {
  return string_shared_field(m_header.field_index(tag));
}

SharedField<int32_t> Variant::shared_field_as_integer(const std::string& tag) const {
  return shared_field_as<int32_t>(tag);
}

SharedField<float> Variant::shared_field_as_float(const std::string& tag) const {
  return shared_field_as<float>(tag);
}

SharedField<string> Variant::shared_field_as_string(const std::string& tag) const {
  return shared_field_as<string>(tag);
}

SharedField<int32_t> Variant::integer_shared_field(const int32_t index) const {
  if (check_field(BCF_HL_INFO, BCF_HT_INT, index))
    return shared_field_as<int32_t>(index); 
  return SharedField<int32_t>{};
}

SharedField<float> Variant::float_shared_field(const int32_t index) const {
  if (check_field(BCF_HL_INFO, BCF_HT_REAL, index))
    return shared_field_as<float>(index); 
  return SharedField<float>{};
}

SharedField<string> Variant::string_shared_field(const int32_t index) const {
  if (check_field(BCF_HL_INFO, BCF_HT_STR, index))
    return shared_field_as<string>(index); 
  return SharedField<string>{};
}

SharedField<int32_t> Variant::shared_field_as_integer(const int32_t index) const {
  return shared_field_as<int32_t>(index);
}

SharedField<float> Variant::shared_field_as_float(const int32_t index) const {
  return shared_field_as<float>(index);
}

SharedField<string> Variant::shared_field_as_string(const int32_t index) const {
  return shared_field_as<string>(index);
}

IndividualField<Genotype> Variant::genotypes() const {
  // bcf_get_fmt() will unpack the record if necessary
  const auto fmt = bcf_get_fmt(m_header.m_header.get(), m_body.get(), "GT");
  if (fmt == nullptr) ///< if the variant is missing or the GT tag is missing, return an empty IndividualField
    return IndividualField<Genotype>{};
  return IndividualField<Genotype>{m_body, fmt};
}

}
