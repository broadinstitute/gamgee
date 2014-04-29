#include "hts_memory.h"

using namespace std;

namespace gamgee {

  /**
   * @brief wraps a pre-allocated bam1_t in a shared_ptr with correct deleter
   */
  shared_ptr<bam1_t> make_shared_bam(bam1_t* bam_ptr) {
    return shared_ptr<bam1_t>(bam_ptr, BamDeleter());
  }

  /**
   * @brief wraps a pre-allocated bam_hdr_t in a shared_ptr with correct deleter
   */
  shared_ptr<bam_hdr_t> make_shared_bam_header(bam_hdr_t* bam_header_ptr) {
    return shared_ptr<bam_hdr_t>(bam_header_ptr, HeaderDeleter());
  }

  /**
   * @brief creates a deep copy of an existing bam1_t
   */
  bam1_t* bam_deep_copy(bam1_t* original) {
    return bam_dup1(original);
  }

  /**
   * @brief creates a deep copy of an existing bam_hdr_t
   */
  bam_hdr_t* bam_header_deep_copy(bam_hdr_t* original) {
    return bam_hdr_dup(original);
  }

}
