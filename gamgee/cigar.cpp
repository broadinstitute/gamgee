#include "cigar.h"
#include "utils/hts_memory.h"

#include <string>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace gamgee {

/**
  * @warning the order of these operators are matching the defines in htslib, do not change!
  */
const char Cigar::cigar_ops_as_chars[] = { 'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B' };

/**
  * @brief creates a Cigar object that points to htslib memory already allocated
  *
  * @note the resulting Cigar object shares ownership of the pre-allocated memory via
  *       shared_ptr reference counting
  */
Cigar::Cigar(const std::shared_ptr<bam1_t>& sam_record) :
  m_sam_record { sam_record },
  m_cigar { bam_get_cigar(sam_record.get()) },
  m_num_cigar_elements { (sam_record.get())->core.n_cigar }
{}

/**
  * @brief creates a deep copy of a Cigar object
  *
  * @note the copy will have exclusive ownership over the newly-allocated htslib memory
  */
Cigar::Cigar(const Cigar& other) :
  m_sam_record { utils::make_shared_sam(utils::sam_deep_copy(other.m_sam_record.get())) },
  m_cigar { bam_get_cigar(m_sam_record.get()) },
  m_num_cigar_elements { other.m_num_cigar_elements }
{}

/**
  * @brief moves a Cigar object, transferring ownership of the underlying htslib memory
  */
Cigar::Cigar(Cigar&& other) noexcept :
  m_sam_record { move(other.m_sam_record) },
  m_cigar { other.m_cigar },
  m_num_cigar_elements { other.m_num_cigar_elements }
{
  other.m_cigar = nullptr;
}

/**
  * @brief creates a deep copy of a Cigar object
  *
  * @note the copy will have exclusive ownership over the newly-allocated htslib memory
  */
Cigar& Cigar::operator=(const Cigar& other) {
  if (&other == this) ///< check for self assignment
    return *this; 
  m_sam_record = utils::make_shared_sam(utils::sam_deep_copy(other.m_sam_record.get())); ///< shared_ptr assignment will take care of deallocating old sam record if necessary
  m_cigar = bam_get_cigar(m_sam_record.get());
  m_num_cigar_elements = other.m_num_cigar_elements;
  return *this;
}

/**
  * @brief moves a Cigar object, transferring ownership of the underlying htslib memory
  */
Cigar& Cigar::operator=(Cigar&& other) noexcept {
  if ( &other == this )  
    return *this;
  m_sam_record = move(other.m_sam_record);
  m_cigar = other.m_cigar;
  other.m_cigar = nullptr;
  m_num_cigar_elements = other.m_num_cigar_elements;
  return *this;
}

/**
  * @brief access an individual cigar element by index
  *
  * @return cigar element at the specified index as an encoded uint32_t. use cigar_op() and
  *         cigar_oplen() to unpack the cigar operator and length
  */
CigarElement Cigar::operator[](const uint32_t index) const {
  if ( index >= m_num_cigar_elements )
    throw out_of_range(string("Index ") + std::to_string(index) + " out of range in Cigar::operator[]");
  return m_cigar[index];
}

/**
  * @brief access and/or modify an individual cigar element by index
  *
  * @return cigar element at the specified index as an encoded uint32_t. use cigar_op() and
  *         cigar_oplen() to unpack the cigar operator and length
  */
CigarElement& Cigar::operator[](const uint32_t index) {
  if ( index >= m_num_cigar_elements )
    throw out_of_range(string("Index ") + std::to_string(index) + " out of range in Cigar::operator[]");
  return m_cigar[index];
}

bool Cigar::operator==(const Cigar& other) const {
  if ( m_num_cigar_elements != other.m_num_cigar_elements )
    return false;

  for ( auto i = 0u; i < m_num_cigar_elements; ++i ) {
    if ( m_cigar[i] != other.m_cigar[i] )
      return false;
  }

  return true;
}

bool Cigar::operator!=(const Cigar& other) const {
  return !(*this == other);
}


/**
  * @brief returns a string representation of this cigar
  */
string Cigar::to_string() const {
  stringstream stream;

  for ( uint32_t i = 0; i < m_num_cigar_elements; ++i ) 
    stream << cigar_oplen(m_cigar[i]) << cigar_ops_as_chars[static_cast<int>(cigar_op(m_cigar[i]))];
  return stream.str();
}


}
