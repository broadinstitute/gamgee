#ifndef __gamgee__base_quals__
#define __gamgee__base_quals__

#include "htslib/sam.h"

#include <memory>

namespace gamgee {

/**
 * @brief Utility class to handle the memory management of the sam record object for a read base qualities
 */
class BaseQuals {
 public:
  explicit BaseQuals(const std::shared_ptr<bam1_t>& sam_record);
  BaseQuals(const BaseQuals& other);
  BaseQuals(BaseQuals&& other) noexcept;
  BaseQuals& operator=(const BaseQuals& other);
  BaseQuals& operator=(BaseQuals&& other) noexcept;
  ~BaseQuals() = default; ///< Default destruction is sufficient, since our shared_ptr will handle deallocation

  uint8_t operator[](const uint32_t index) const; ///< use freely as you would an array. @note currently implemented as read only
  uint32_t size() const { return m_num_quals; }   ///< number of base qualities in the container

 private:
  std::shared_ptr<bam1_t> m_sam_record; ///< sam record containing our base qualities, potentially co-owned by multiple other objects
  uint8_t* m_quals;                     ///< Pointer to the start of the base qualities in m_sam_record, cached for efficiency
  uint32_t m_num_quals;                 ///< Number of quality scores in our sam record
};

}

#endif /* __gamgee__base_quals__ */
