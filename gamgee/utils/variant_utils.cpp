#include "variant_utils.h"

#include "htslib/vcf.h"

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
}

}
