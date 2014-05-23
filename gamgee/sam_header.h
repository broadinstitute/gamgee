#ifndef __gamgee__sam_header__
#define __gamgee__sam_header__

#include "htslib/sam.h"

#include <memory>

namespace gamgee {

/**
 * @brief Utility class to hold the header of a sam file
 */
class SamHeader {
 public:
  explicit SamHeader() = default;                               ///< @brief initializes a null SamHeader @warning if you need to create a SamHeader from scratch, use the builder instead
  explicit SamHeader(const std::shared_ptr<bam_hdr_t>& header); ///< @brief creates a SamHeader given htslib object. @note used by all iterators
  SamHeader(const SamHeader& other);                            ///< @brief makes a deep copy of a SamHeader. Shared pointers maintain state to all other associated objects correctly.
  SamHeader(SamHeader&& other) noexcept;                        ///< @brief moves SamHeader accordingly. Shared pointers maintain state to all other associated objects correctly.
  SamHeader& operator=(const SamHeader& other);                 ///< @brief deep copy assignment of a SamHeader. Shared pointers maintain state to all other associated objects correctly.
  SamHeader& operator=(SamHeader&& other) noexcept;             ///< @brief move assignment of a SamHeader. Shared pointers maintain state to all other associated objects correctly.

 private:
  std::shared_ptr<bam_hdr_t> m_header;

  friend class SamWriter;
};

}
#endif 
