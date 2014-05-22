#ifndef __gamgee__sam__
#define __gamgee__sam__

#include "sam_header.h"
#include "sam_body.h"

#include "htslib/sam.h"

#include <string>
#include <memory>

namespace gamgee {

/**
 * @brief Utility class to manipulate a Sam record.
 */
class Sam : public SamBody {
 public:
  explicit Sam() = default;
  explicit Sam(const std::shared_ptr<bam1_t>& body, const std::shared_ptr<bam_hdr_t>& header) noexcept : SamBody{body}, m_header{header} {}
  SamHeader header() {
    return SamHeader{m_header}; ///< TODO: return a reference to a singleton SamHeader object instead of constructing a new one each time, since we're already sharing the underlying htslib memory
  }

 private:
  std::shared_ptr<bam_hdr_t> m_header;
};

}  // end of namespace

#endif /* defined(__gamgee__sam__) */
