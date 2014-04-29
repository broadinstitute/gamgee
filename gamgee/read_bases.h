#ifndef __gamgee__read_bases__
#define __gamgee__read_bases__

#include "htslib/sam.h"

#include <memory>
#include <map>

namespace gamgee {

// Enum values used here correspond to the 4-bit base encodings in htslib so that we can
// cast directly to Base
enum class Base { A = 1, C = 2, G = 4, T = 8, N = 15 };

class ReadBases {
public:
  explicit ReadBases(const std::shared_ptr<bam1_t>& sam_record);
  ReadBases(const ReadBases& other);
  ReadBases(ReadBases&& other);
  ReadBases& operator=(const ReadBases& other);
  ReadBases& operator=(ReadBases&& other);

  // Default destruction is sufficient, since our shared_ptr will handle deallocation
  ~ReadBases() = default;

  Base operator[](const uint32_t index) const;
  uint32_t size() const { return m_num_bases; };
  std::string to_string() const;

private:
  // Sam record containing our bases, potentially co-owned by multiple other objects
  std::shared_ptr<bam1_t> m_sam_record;

  // Pointer to the start of the bases in m_sam_record, cached for efficiency
  uint8_t* m_bases;

  // Number of bases in our sam record
  uint32_t m_num_bases;

  static const std::map<Base, const char*> base_to_string_map;
};

}

#endif /* __gamgee__read_bases__ */
