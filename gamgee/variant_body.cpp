#include "variant_body.h"
#include "utils/hts_memory.h"

#include "htslib/vcf.h"

using namespace std;

namespace gamgee {

/**
 * @brief creates an empty variant record and allocates new htslib memory for all the fields
 */
VariantBody::VariantBody() :
  m_body {utils::make_shared_variant(bcf_init1())}
{}

/**
 * @brief creates a variant record that points to htslib memory already allocated
 *
 * @note the resulting VariantBody shares ownership of the pre-allocated memory via shared_ptr
 *       reference counting
 */
VariantBody::VariantBody(const std::shared_ptr<bcf1_t>& body) :
  m_body { body }
{}

/**
 * @brief creates a deep copy of a variant record
 */
VariantBody::VariantBody(const VariantBody& other) :
  m_body { utils::make_shared_variant(utils::variant_deep_copy(other.m_body.get())) }
{}

/**
 * @brief moves a variant record, transferring ownership of the underlying htslib memory
 */
VariantBody::VariantBody(VariantBody&& other) noexcept :
  m_body { move(other.m_body) }
{}

/**
 * @brief creates a deep copy of a variant record
 * @param other the VariantBody to be copied
 * @note the parameter is passed by copy so we can maintain the strong exception safety guarantee
 */
VariantBody& VariantBody::operator=(const VariantBody& other) {
  if ( &other == this )  
    return *this;
  m_body = utils::make_shared_variant(utils::variant_deep_copy(other.m_body.get())); ///< shared_ptr assignment will take care of deallocating old sam record if necessary
  return *this;
}

/**
 * @brief moves a variant record, transferring ownership of the underlying htslib memory
 */
VariantBody& VariantBody::operator=(VariantBody&& other) noexcept {
  if ( &other == this )  
    return *this;
  m_body = move(other.m_body);
  return *this;
}


}

