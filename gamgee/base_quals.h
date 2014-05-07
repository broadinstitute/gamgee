#ifndef __gamgee__base_quals__
#define __gamgee__base_quals__

#include "htslib/sam.h"

#include <memory>

namespace gamgee {

class BaseQuals {
public:
  explicit BaseQuals(const std::shared_ptr<bam1_t>& sam_record);
  BaseQuals(const BaseQuals& other);
  BaseQuals(BaseQuals&& other);
  BaseQuals& operator=(const BaseQuals& other);
  BaseQuals& operator=(BaseQuals&& other);

  // Default destruction is sufficient, since our shared_ptr will handle deallocation
  ~BaseQuals() = default;

  uint8_t operator[](const uint32_t index) const;
  uint32_t size() const { return m_num_quals; }

private:
  // Sam record containing our base qualities, potentially co-owned by multiple other objects
  std::shared_ptr<bam1_t> m_sam_record;

  // Pointer to the start of the base qualities in m_sam_record, cached for efficiency
  uint8_t* m_quals;

  // Number of quality scores in our sam record
  uint32_t m_num_quals;
};

}

#endif /* __gamgee__base_quals__ */
