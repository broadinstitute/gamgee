#include "base_quals.h"

#include "../utils/hts_memory.h"
#include "../utils/utils.h"

#include <string>
#include <stdexcept>
#include <sstream>

using namespace std;

namespace gamgee {

/**
  * @brief creates a BaseQuals object that points to htslib memory already allocated
  *
  * @note the resulting BaseQuals object shares ownership of the pre-allocated memory via
  *       shared_ptr reference counting
  */
BaseQuals::BaseQuals(const std::shared_ptr<bam1_t>& sam_record) :
  m_sam_record { sam_record },
  m_quals { bam_get_qual(sam_record.get()) },
  m_num_quals { uint32_t((sam_record.get())->core.l_qseq) }
{}

/**
  * @brief creates a deep copy of a BaseQuals object
  *
  * @note the copy will have exclusive ownership over the newly-allocated htslib memory
  */
BaseQuals::BaseQuals(const BaseQuals& other) :
  m_sam_record { utils::make_shared_sam(utils::sam_deep_copy(other.m_sam_record.get())) },
  m_quals { bam_get_qual(m_sam_record.get()) },
  m_num_quals { other.m_num_quals }
{}

/**
  * @brief creates a deep copy of a BaseQuals object
  *
  * @note the copy will have exclusive ownership over the newly-allocated htslib memory
  */
BaseQuals& BaseQuals::operator=(const BaseQuals& other) {
  if ( &other == this )  
    return *this;
  m_sam_record = utils::make_shared_sam(utils::sam_deep_copy(other.m_sam_record.get())); ///< shared_ptr assignment will take care of deallocating old sam record if necessary
  m_quals = bam_get_qual(m_sam_record.get()); 
  m_num_quals = other.m_num_quals;
  return *this;
}

/**
  * @brief access an individual base quality by index
  *
  * @return base quality at the specified index as an unsigned byte
  */
uint8_t BaseQuals::operator[](const uint32_t index) const {
  utils::check_max_boundary(index, m_num_quals);
  return m_quals[index];
}

/**
 * @brief access and/or modify an individual base quality by index
 *
 * @return base quality at the specified index as an unsigned byte
 */
uint8_t& BaseQuals::operator[](const uint32_t index) {
  utils::check_max_boundary(index, m_num_quals);
  return m_quals[index];
}

/**
 * @brief check whether this object contains the same base qualities as another BaseQuals object
 */
bool BaseQuals::operator==(const BaseQuals& other) const {
  if ( m_num_quals != other.m_num_quals )
    return false;

  for ( auto i = 0u; i < m_num_quals; ++i ) {
    if ( m_quals[i] != other.m_quals[i] )
      return false;
  }

  return true;
}

/**
 * @brief check whether this object does not contain the same base qualities as another BaseQuals object
 */
bool BaseQuals::operator!=(const BaseQuals& other) const {
  return !(*this == other);
}

/**
 * @brief produce a string representation of the base qualities in this object
 */
std::string BaseQuals::to_string() const {
  stringstream stream;

  for ( auto i = 0u; i < m_num_quals; ++i ) {
    stream << int(m_quals[i]);
    if ( i < m_num_quals - 1 )
      stream << " ";
  }
  return stream.str();
}

} // end of namespace gamgee
