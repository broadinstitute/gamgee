#include "variant_utils.h"

#include "../exceptions.h"
#include "hts_memory.h"

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
    std::for_each(samples.begin(), samples.end(), [&sample_list](const auto& s) { sample_list += s + ","; });
    sample_list.erase(sample_list.size() - 1);
    bcf_hdr_set_samples(hdr_ptr, sample_list.c_str(), false);
  }

  // NOTE: must NOT call bcf_hdr_sync() here, since htslib calls it for us in bcf_hdr_set_samples()
}

void merge_variant_headers(const std::shared_ptr<bcf_hdr_t>& dest_hdr_ptr, const std::shared_ptr<bcf_hdr_t>& src_hdr_ptr) {
  auto success = bcf_hdr_combine(dest_hdr_ptr.get(), src_hdr_ptr.get());
  if (success != 0)
    throw HtslibException(success);

  // TODO: there is probably a more efficient way
  for (auto sample_counter = 0; sample_counter < bcf_hdr_nsamples(src_hdr_ptr.get()); ++sample_counter) {
    // don't check for error code because the only "error" is ignoring a duplicate sample, not an error for us
    bcf_hdr_add_sample(dest_hdr_ptr.get(), src_hdr_ptr->samples[sample_counter]);
  }

  // vcf.h    "After all samples have been added, NULL must be passed to update internal header structures."
  bcf_hdr_add_sample(dest_hdr_ptr.get(), nullptr);
  bcf_hdr_sync(dest_hdr_ptr.get());
}

}
