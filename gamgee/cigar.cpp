#include "cigar.h"
#include "hts_memory.h"

#include <string>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace gamgee {

  /**
   * @brief creates a Cigar object that points to htslib memory already allocated
   *
   * @note the resulting Cigar object shares ownership of the pre-allocated memory via
   *       shared_ptr reference counting
   */
  Cigar::Cigar(const shared_ptr<bam1_t>& sam_record) :
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
    m_sam_record { make_shared_bam(bam_deep_copy(other.m_sam_record.get())) },
    m_cigar { bam_get_cigar(m_sam_record.get()) },
    m_num_cigar_elements { other.m_num_cigar_elements }
  {}

  /**
   * @brief moves a Cigar object, transferring ownership of the underlying htslib memory
   */
  Cigar::Cigar(Cigar&& other) :
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
    if ( &other != this ) {
      // shared_ptr assignment will take care of deallocating old sam record if necessary
      m_sam_record = make_shared_bam(bam_deep_copy(other.m_sam_record.get()));
      m_cigar = bam_get_cigar(m_sam_record.get());
      m_num_cigar_elements = other.m_num_cigar_elements;
    }

    return *this;
  }

  /**
   * @brief moves a Cigar object, transferring ownership of the underlying htslib memory
   */
  Cigar& Cigar::operator=(Cigar&& other) {
    if ( &other != this ) {
      m_sam_record = move(other.m_sam_record);
      m_cigar = other.m_cigar;
      other.m_cigar = nullptr;
      m_num_cigar_elements = other.m_num_cigar_elements;
    }

    return *this;
  }

  /**
   * @brief access an individual cigar element by index
   *
   * @return cigar element at the specified index as an encoded uint32_t. use cigar_op() and
   *         cigar_oplen() to unpack the cigar operator and length
   */
  uint32_t Cigar::operator[](const uint32_t index) const {
    if ( index >= m_num_cigar_elements )
      throw out_of_range(string("Index ") + std::to_string(index) + " out of range in Cigar::operator[]");

    return m_cigar[index];
  }

  // Lookup table to convert CigarOperator enum values to chars. An unfortunate necessity due
  // to lack of features in C++ enum "classes"
  const char Cigar::cigar_ops_as_chars[] = { 'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B' };

  /**
   * @brief returns a string representation of this cigar
   */
  string Cigar::to_string() const {
    stringstream stream;

    for ( uint32_t i = 0; i < m_num_cigar_elements; ++i ) {
      stream << cigar_oplen(m_cigar[i]) << cigar_ops_as_chars[static_cast<int>(cigar_op(m_cigar[i]))];
    }

    return stream.str();
  }
}
