#include "sam_header.h"
#include "hts_memory.h"

using namespace std;

namespace gamgee {

  /**
   * @brief creates an empty SamHeader, allocating new htslib memory
   */
  SamHeader::SamHeader() :
    m_header { make_shared_bam_header(bam_hdr_init()) }
  {}

  /**
   * @brief creates a SamHeader object that points to htslib memory already allocated
   *
   * @note the resulting SamHeader object shares ownership of the pre-allocated memory via
   *       shared_ptr reference counting
   */
  SamHeader::SamHeader(const shared_ptr<bam_hdr_t>& header) :
    m_header { header }
  {}

  /**
   * @brief creates a deep copy of a SamHeader object
   *
   * @note the copy will have exclusive ownership over the newly-allocated htslib memory
   */
  SamHeader::SamHeader(const SamHeader& other) :
    m_header { make_shared_bam_header(bam_header_deep_copy(other.m_header.get())) }
  {}

  /**
   * @brief moves a SamHeader object, transferring ownership of the underlying htslib memory
   */
  SamHeader::SamHeader(SamHeader&& other) :
    m_header { move(other.m_header) }
  {}

  /**
   * @brief creates a deep copy of a SamHeader object
   *
   * @note the copy will have exclusive ownership over the newly-allocated htslib memory
   */
  SamHeader& SamHeader::operator=(const SamHeader& other) {
    if ( &other != this ) {
      // shared_ptr assignment will take care of deallocating old sam header if necessary
      m_header = make_shared_bam_header(bam_header_deep_copy(other.m_header.get()));
    }

    return *this;
  }

  SamHeader& SamHeader::operator=(SamHeader&& other) {
    if ( &other != this ) {
      // shared_ptr assignment will take care of deallocating old sam header if necessary
      m_header = move(other.m_header);
    }

    return *this;
  }

}
