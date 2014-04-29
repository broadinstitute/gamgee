#include "read_bases.h"
#include "hts_memory.h"

#include <string>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace gamgee {

  /**
   * @brief creates a ReadBases object that points to htslib memory already allocated
   *
   * @note the resulting ReadBases object shares ownership of the pre-allocated memory via
   *       shared_ptr reference counting
   */
  ReadBases::ReadBases(const shared_ptr<bam1_t>& sam_record) :
    m_sam_record { sam_record },
    m_bases { bam_get_seq(sam_record.get()) },
    m_num_bases { uint32_t((sam_record.get())->core.l_qseq) }
  {}

  /**
   * @brief creates a deep copy of a ReadBases object
   *
   * @note the copy will have exclusive ownership over the newly-allocated htslib memory
   */
  ReadBases::ReadBases(const ReadBases& other) :
    m_sam_record { make_shared_bam(bam_deep_copy(other.m_sam_record.get())) },
    m_bases { bam_get_seq(m_sam_record.get()) },
    m_num_bases { other.m_num_bases }
  {}

  /**
   * @brief moves a ReadBases object, transferring ownership of the underlying htslib memory
   */
  ReadBases::ReadBases(ReadBases&& other) :
    m_sam_record { move(other.m_sam_record) },
    m_bases { other.m_bases },
    m_num_bases { other.m_num_bases }
  {
    other.m_bases = nullptr;
  }

  /**
   * @brief creates a deep copy of a ReadBases object
   *
   * @note the copy will have exclusive ownership over the newly-allocated htslib memory
   */
  ReadBases& ReadBases::operator=(const ReadBases& other) {
    if ( &other != this ) {
      // shared_ptr assignment will take care of deallocating old sam record if necessary
      m_sam_record = make_shared_bam(bam_deep_copy(other.m_sam_record.get()));
      m_bases = bam_get_seq(m_sam_record.get());
      m_num_bases = other.m_num_bases;
    }

    return *this;
  }

  /**
   * @brief moves a ReadBases object, transferring ownership of the underlying htslib memory
   */
  ReadBases& ReadBases::operator=(ReadBases&& other) {
    if ( &other != this ) {
      // shared_ptr assignment will take care of deallocating old sam record if necessary
      m_sam_record = move(other.m_sam_record);
      m_bases = other.m_bases;
      other.m_bases = nullptr;
      m_num_bases = other.m_num_bases;
    }

    return *this;
  }

  /**
   * @brief access an individual base by index
   *
   * @return base at the specified index as an enumerated value
   */
  Base ReadBases::operator[](const uint32_t index) const {
    if ( index >= m_num_bases )
      throw out_of_range(string("Index ") + std::to_string(index) + " out of range in ReadBases::operator[]");

    return static_cast<Base>(bam_seqi(m_bases, index));
  }

  // Lookup table to convert Base enum values to chars. An unfortunate necessity due
  // to lack of features in C++ enum "classes"
  const map<Base, const char*> ReadBases::base_to_string_map = { {Base::A, "A"}, {Base::C, "C"}, {Base::G, "G"}, {Base::T, "T"}, {Base::N, "N"} };

  /**
   * @brief returns a string representation of the bases in this read
   */
  string ReadBases::to_string() const {
    stringstream ss;

    for ( uint32_t i = 0; i < m_num_bases; i++ ) {
      ss << base_to_string_map.at((*this)[i]);
    }

    return ss.str();
  }
}
