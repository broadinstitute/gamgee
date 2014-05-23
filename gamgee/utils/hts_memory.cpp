#include "hts_memory.h"

#include <memory>

using namespace std;

namespace gamgee {
namespace utils {

/**
  * @brief wraps a pre-allocated bam1_t in a shared_ptr with correct deleter
  * @param sam_ptr an htslib raw bam pointer
  */
shared_ptr<bam1_t> make_shared_sam(bam1_t* sam_ptr) {
  return shared_ptr<bam1_t>(sam_ptr, SamBodyDeleter());
}

/**
  * @brief wraps a pre-allocated bam_hdr_t in a shared_ptr with correct deleter
  * @param sam_header_ptr an htslib raw bam header pointer
  */
shared_ptr<bam_hdr_t> make_shared_sam_header(bam_hdr_t* sam_header_ptr) {
  return shared_ptr<bam_hdr_t>(sam_header_ptr, SamHeaderDeleter());
}

/**
  * @brief wraps a pre-allocated bcf1_t in a shared_ptr with correct deleter
  * @param bcf_ptr an htslib raw vcf pointer
  */
shared_ptr<bcf1_t> make_shared_variant(bcf1_t* bcf_ptr) {
  return shared_ptr<bcf1_t>(bcf_ptr, VariantBodyDeleter());
}

/**
  * @brief wraps a pre-allocated bcf_hdr_t in a shared_ptr with correct deleter
  * @param bcf_hdr_ptr an htslib raw variant header pointer
  */
shared_ptr<bcf_hdr_t> make_shared_variant_header(bcf_hdr_t* bcf_hdr_ptr) {
  return shared_ptr<bcf_hdr_t>(bcf_hdr_ptr, VariantHeaderDeleter());
}

/**
  * @brief creates a deep copy of an existing bam1_t
  * @param original an htslib raw bam pointer
  */
bam1_t* sam_deep_copy(bam1_t* original) {
  return bam_dup1(original);
}

/**
  * @brief creates a deep copy of an existing bam_hdr_t
  * @param original an htslib raw bam header pointer
  */
bam_hdr_t* sam_header_deep_copy(bam_hdr_t* original) {
  return bam_hdr_dup(original);
}

/**
 * @brief creates a deep copy of an existing bcf1_t
 * @param original an htslib raw bcf pointer
 */
bcf1_t* variant_deep_copy(bcf1_t* original) {
  return bcf_dup(original);
}

/**
  * @brief creates a deep copy of an existing bcf_hdr_t
  * @param original an htslib raw bcf header pointer
  */
bcf_hdr_t* variant_header_deep_copy(bcf_hdr_t* original) {
  return bcf_hdr_dup(original);
}



}
}
