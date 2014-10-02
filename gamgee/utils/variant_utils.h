#ifndef gamgee__variant_utils__guard
#define gamgee__variant_utils__guard

#include "hts_memory.h"

#include "htslib/vcf.h"

#include <string>
#include <vector>

namespace gamgee {

/**
 * @brief allows the caller to include only selected samples in a Variant Reader. To create a sites only file,
 * simply pass an empty vector of samples.
 *
 * @param samples the list of samples you want included/excluded from your iteration
 * @param whether you want these samples to be included or excluded from your iteration.
 */
void subset_variant_samples(bcf_hdr_t* hdr_ptr, const std::vector<std::string>& samples, const bool include);

enum class AlleleType { REFERENCE, SNP, INSERTION, DELETION };

using AlleleMask = std::vector<AlleleType>;

/**
 * @brief merges a variant header into another
 *
 * @param dest_hdr_ptr a shared pointer to a bcf_hdr_t containing the header to be merged into
 * @param src_hdr_ptr a shared pointer to a bcf_hdr_t containing the header to be merged from
 */
void merge_variant_headers(const std::shared_ptr<bcf_hdr_t>& dest_hdr_ptr, const std::shared_ptr<bcf_hdr_t>& src_hdr_ptr);
}

#endif /* gamgee__variant_utils__guard */
