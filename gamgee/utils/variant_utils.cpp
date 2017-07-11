#include "variant_utils.h"

#include "hts_memory.h"

#include "../exceptions.h"

#include "htslib/vcf.h"

#include <algorithm>
#include <string>
#include <vector>

namespace gamgee {

void subset_variant_samples(bcf_hdr_t* hdr_ptr, const std::vector<std::string>& samples, const bool include) {
  if (samples.empty() && include) // exclude all samples
    bcf_hdr_set_samples(hdr_ptr, NULL, false);

  else if (samples.empty() && !include) // keep all samples
    bcf_hdr_set_samples(hdr_ptr, "-", false);

  else { // select some samples
    auto sample_list = include ? std::string{} : std::string{"^"};
    std::for_each(samples.begin(), samples.end(), [&sample_list](const std::string& s) { sample_list += s + ","; });
    sample_list.erase(sample_list.size() - 1);
    bcf_hdr_set_samples(hdr_ptr, sample_list.c_str(), false);
  }

  // NOTE: must NOT call bcf_hdr_sync() here, since htslib calls it for us in bcf_hdr_set_samples()
}

void merge_variant_headers(const std::shared_ptr<bcf_hdr_t>& dest_hdr_ptr, const std::shared_ptr<bcf_hdr_t>& src_hdr_ptr) {
  bcf_hdr_merge(dest_hdr_ptr.get(), src_hdr_ptr.get());
  // TODO: there is probably a more efficient way
  for (auto sample_counter = 0; sample_counter < bcf_hdr_nsamples(src_hdr_ptr.get()); ++sample_counter) {
    // Note: bcf_hdr_add_sample aborts if there are duplicated sample names.
    bcf_hdr_add_sample(dest_hdr_ptr.get(), src_hdr_ptr->samples[sample_counter]);
  }

  bcf_hdr_sync(dest_hdr_ptr.get());
}

}
