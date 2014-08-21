#ifndef gamgee__hts_memory__guard
#define gamgee__hts_memory__guard

#include "htslib/sam.h"
#include "htslib/vcf.h"

#include <memory>
#include <vector>
#include <string>

namespace gamgee {
namespace utils { 

/** 
 * @brief a functor object to delete an htsFile pointer
 */
struct HtsFileDeleter {
  void operator()(htsFile* p) const { hts_close(p); }
};

/** 
 * @brief a functor object to delete an hts file index pointer
 */
struct HtsIndexDeleter {
  void operator()(hts_idx_t* p) const { hts_idx_destroy(p); }
};

/**
 * @brief a functor object to delete an hts file iterator pointer
 */
struct HtsIteratorDeleter {
  void operator()(hts_itr_t* p) const { hts_itr_destroy(p); }
};

/**
 * @brief a functor object to delete a bam1_t pointer 
 * 
 * used by the unique_ptr queue used in the paired reader
 */
struct SamBodyDeleter {
  void operator()(bam1_t* p) const { bam_destroy1(p); }
};

/** 
 * @brief a functor object to delete a bam_hdr_t pointer 
 */
struct SamHeaderDeleter {
  void operator()(bam_hdr_t* p) const { bam_hdr_destroy(p); }
};

/** 
 * @brief a functor object to delete a bcf_hdr_t pointer
 */
struct VariantHeaderDeleter {
  void operator()(bcf_hdr_t* p) const { bcf_hdr_destroy(p); }
};

/** 
 * @brief a functor object to delete a bcf1_t pointer
 */
struct VariantBodyDeleter {
  void operator()(bcf1_t* p) const { bcf_destroy1(p); }
};

std::shared_ptr<htsFile> make_shared_hts_file(htsFile* hts_file_ptr);
std::shared_ptr<hts_idx_t> make_shared_hts_index(hts_idx_t* hts_index_ptr);
std::shared_ptr<hts_itr_t> make_shared_hts_itr(hts_itr_t* hts_itr_ptr);
std::shared_ptr<bam1_t> make_shared_sam(bam1_t* sam_ptr);
std::shared_ptr<bam_hdr_t> make_shared_sam_header(bam_hdr_t* sam_header_ptr);
std::shared_ptr<bcf1_t> make_shared_variant(bcf1_t* bcf_ptr);
std::shared_ptr<bcf_hdr_t> make_shared_variant_header(bcf_hdr_t* bcf_hdr_ptr);

bam1_t* sam_deep_copy(bam1_t* original);
bam_hdr_t* sam_header_deep_copy(bam_hdr_t* original);
bcf1_t* variant_deep_copy(bcf1_t* original); 
bcf_hdr_t* variant_header_deep_copy(bcf_hdr_t* original);

bam1_t* sam_shallow_copy(bam1_t* original);

/**
 * @brief helper function to translate an index into a string in the filter list 
 * @param header a VariantHeader htslib pointer
 * @param body a Variant htslib pointer
 * @param index the index of the filter you want access to
 */
std::string htslib_filter_name(bcf_hdr_t* header, bcf1_t* body, int index);

}
}

#endif // gamgee__hts_memory__guard
