#include "variant.h"
#include "utils/hts_memory.h"

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

vector<uint8_t> Variant::genotype_quals() const {
  const auto fmt = bcf_get_fmt(m_header.get(), m_body.get(), "GQ");     
  if (fmt == nullptr) ///< if the format is missing or the GQ tag is missing, return an empty vector
    return vector<uint8_t>{};
  return vector<uint8_t>(fmt->p, fmt->p + fmt->p_len); 
}

}

