#include "hts_memory.h"

#include <memory>

using namespace std;

namespace gamgee {
namespace utils {

/**
  * @brief wraps a pre-allocated htsFile in a shared_ptr with correct deleter
  * @param hts_file_ptr an htslib raw file pointer
  */
shared_ptr<htsFile> make_shared_hts_file(htsFile* hts_file_ptr) {
  return shared_ptr<htsFile>(hts_file_ptr, HtsFileDeleter());
}

/**
  * @brief wraps a pre-allocated hts_idx_t in a shared_ptr with correct deleter
  * @param hts_index_ptr an htslib raw file index pointer
  */
shared_ptr<hts_idx_t> make_shared_hts_index(hts_idx_t* hts_index_ptr) {
  return shared_ptr<hts_idx_t>(hts_index_ptr, HtsIndexDeleter());
}

/**
  * @brief wraps a pre-allocated hts_itr_t in a shared_ptr with correct deleter
  * @param hts_itr_ptr an htslib raw file iterator pointer
  */
shared_ptr<hts_itr_t> make_shared_hts_itr(hts_itr_t* hts_itr_ptr) {
  return shared_ptr<hts_itr_t>(hts_itr_ptr, HtsIteratorDeleter());
}

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

/**
 * @brief creates a shallow copy of an existing bam1_t: copies
 *        core fields but not the data buffer or fields related
 *        to the size of the data buffer
 */
bam1_t* sam_shallow_copy(bam1_t* original) {
  bam1_t* new_read = bam_init1();

  // Copy full struct contents, then zero out fields related to the data
  *new_read = *original;

  new_read->data = nullptr;
  new_read->l_data = 0;
  new_read->m_data = 0;
  new_read->core.l_qname = 0;
  new_read->core.l_qseq = 0;
  new_read->core.n_cigar = 0;

  return new_read;
}

std::string htslib_filter_name(bcf_hdr_t* header, bcf1_t* body, int index) { 
  return std::string{bcf_hdr_int2id(header, BCF_DT_ID, body->d.flt[index])};
} 

}
}
