#ifndef __gamgee__read_bases__
#define __gamgee__read_bases__

#include "htslib/sam.h"

#include <memory>
#include <map>

namespace gamgee {

/**
 * @brief simple enum to hold all valid bases in the SAM format
 * @note enum values used here correspond to the 4-bit base encodings in htslib 
 * so that we can cast directly to Base 
 */
enum class Base { A = 1, C = 2, G = 4, T = 8, N = 15 };

/**
 * @brief Utility class to handle the memory management of the sam record object for read bases
 *
 * This class uses Base to represent the bases A,C,G,T,N so we can get byte by byte correspondence
 * with the underlying compressed memory model. 
 *
 * @note Any functionality lost because of this should be made available by the ReadBases class.
 */
class ReadBases {
public:
  explicit ReadBases(const std::shared_ptr<bam1_t>& sam_record);
  ReadBases(const ReadBases& other);
  ReadBases(ReadBases&& other) noexcept;
  ReadBases& operator=(const ReadBases& other);
  ReadBases& operator=(ReadBases&& other) noexcept;
  ~ReadBases() = default; ///< default destruction is sufficient, since our shared_ptr will handle deallocation

  Base operator[](const uint32_t index) const;   ///< use freely as you would an array. @note currently implemented as read only
  void set_base(const uint32_t index, const Base base);  ///< modify a base at a specific index
  uint32_t size() const { return m_num_bases; }; ///< number of base qualities in the container
  bool operator==(const ReadBases& other) const; ///< check for equality with another ReadBases object
  bool operator!=(const ReadBases& other) const; ///< check for inequality with another ReadBases object
  std::string to_string() const;  ///< produce a string representation of the bases in this object

private:
  std::shared_ptr<bam1_t> m_sam_record; ///< sam record containing our bases, potentially co-owned by multiple other objects
  uint8_t* m_bases;                     ///< pointer to the start of the bases in m_sam_record, cached for efficiency
  uint32_t m_num_bases;                 ///< number of bases in our sam record

  static const std::map<Base, const char*> base_to_string_map; ///< @brief simple lookup table to convert Base enum values to chars. 

  friend class SamBuilder; ///< builder needs access to the internals in order to build efficiently
};

}

#endif /* __gamgee__read_bases__ */
