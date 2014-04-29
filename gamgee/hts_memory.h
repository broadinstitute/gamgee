#ifndef __gamgee_hts_memory__
#define __gamgee_hts_memory__

#include "sam.h"

/**
 * @brief private namespace to contain htslib specific wrappers
 */
namespace gamgee {

/** 
 * @brief a functor object to delete a bam1_t pointer 
 * 
 * used by the unique_ptr queue used in the paired reader
 */
struct BamDeleter {
  void operator()(bam1_t* p) const { bam_destroy1(p); }
};

/** 
 * @brief a functor object to delete a bam_hdr_t pointer 
 */
struct HeaderDeleter {
  void operator()(bam_hdr_t* p) const { bam_hdr_destroy(p); }
};

std::shared_ptr<bam1_t> make_shared_bam(bam1_t* bam_ptr);
std::shared_ptr<bam_hdr_t> make_shared_bam_header(bam_hdr_t* bam_header_ptr);

bam1_t* bam_deep_copy(bam1_t* original);
bam_hdr_t* bam_header_deep_copy(bam_hdr_t* original);

}


#endif
