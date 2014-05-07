#ifndef __gamgee__sam_header__
#define __gamgee__sam_header__

#include "htslib/sam.h"

#include <memory>

namespace gamgee {

class SamHeader {
 public:
  explicit SamHeader();
  SamHeader(const std::shared_ptr<bam_hdr_t>& header);
  SamHeader(const SamHeader& other);
  SamHeader(SamHeader&& other);
  SamHeader& operator=(const SamHeader& other);
  SamHeader& operator=(SamHeader&& other);

  // Default destruction is sufficient, since our shared_ptr will handle deallocation
  ~SamHeader() = default;

 private:
  std::shared_ptr<bam_hdr_t> m_header;

  friend class SamWriter;
};

}
#endif 
