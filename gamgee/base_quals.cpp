#include "base_quals.h"
#include "hts_memory.h"

#include <string>
#include <stdexcept>

using namespace std;

namespace gamgee {

  /**
   * @brief creates a BaseQuals object that points to htslib memory already allocated
   *
   * @note the resulting BaseQuals object shares ownership of the pre-allocated memory via
   *       shared_ptr reference counting
   */
  BaseQuals::BaseQuals(const shared_ptr<bam1_t>& sam_record) :
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
    m_sam_record { make_shared_bam(bam_deep_copy(other.m_sam_record.get())) },
    m_quals { bam_get_qual(m_sam_record.get()) },
    m_num_quals { other.m_num_quals }
  {}

  /**
   * @brief moves a BaseQuals object, transferring ownership of the underlying htslib memory
   */
  BaseQuals::BaseQuals(BaseQuals&& other) :
    m_sam_record { move(other.m_sam_record) },
    m_quals { other.m_quals },
    m_num_quals { other.m_num_quals }
  {
    other.m_quals = nullptr;
  }

  /**
   * @brief creates a deep copy of a BaseQuals object
   *
   * @note the copy will have exclusive ownership over the newly-allocated htslib memory
   */
  BaseQuals& BaseQuals::operator=(const BaseQuals& other) {
    if ( &other != this ) {
      // shared_ptr assignment will take care of deallocating old sam record if necessary
      m_sam_record = make_shared_bam(bam_deep_copy(other.m_sam_record.get()));
      m_quals = bam_get_qual(m_sam_record.get());
      m_num_quals = other.m_num_quals;
    }

    return *this;
  }

  /**
   * @brief moves a BaseQuals object, transferring ownership of the underlying htslib memory
   */
  BaseQuals& BaseQuals::operator=(BaseQuals&& other) {
    if ( &other != this ) {
      // shared_ptr assignment will take care of deallocating old sam record if necessary
      m_sam_record = move(other.m_sam_record);
      m_quals = other.m_quals;
      other.m_quals = nullptr;
      m_num_quals = other.m_num_quals;
    }

    return *this;
  }

  /**
   * @brief access an individual base quality by index
   *
   * @return base quality at the specified index as an unsigned byte
   */
  uint8_t BaseQuals::operator[](const uint32_t index) const {
    if ( index >= m_num_quals )
      throw out_of_range(string("Index ") + to_string(index) + " out of range in BaseQuals::operator[]");

    return m_quals[index];
  }
}
