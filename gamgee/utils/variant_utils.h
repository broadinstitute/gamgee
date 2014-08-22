#ifndef gamgee__variant_utils__guard
#define gamgee__variant_utils__guard

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

}

#endif /* gamgee__variant_utils__guard */
