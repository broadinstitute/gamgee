#ifndef __gamgee__cigar__
#define __gamgee__cigar__

#include "htslib/sam.h"

#include <memory>
#include <string>

namespace gamgee {

/**
 * @brief comprehensive list of valid cigar operators
 * @note order of the operators in this enum must match the order in BAM_CIGAR_STR from htslib/sam.h
 */
enum class CigarOperator { M, I, D, N, S, H, P, EQ, X, B };

/**
 * @brief Utility class to manage the memory of the cigar structure
 */
class Cigar {
 public:
  explicit Cigar(const std::shared_ptr<bam1_t>& sam_record);
  Cigar(const Cigar& other);
  Cigar(Cigar&& other) noexcept;
  Cigar& operator=(const Cigar& other);
  Cigar& operator=(Cigar&& other) noexcept;
  ~Cigar() = default; ///< default destruction is sufficient, since our shared_ptr will handle deallocation

  uint32_t operator[](const uint32_t index) const;       ///< use freely as you would an array. @note currently implemented as read only
  uint32_t size() const { return m_num_cigar_elements; } ///< number of base qualities in the container
  std::string to_string() const;

  /**
   * @brief gets the operator of an individual cigar element
   */
  inline static CigarOperator cigar_op(const uint32_t cigar_element) {
    return static_cast<CigarOperator>(bam_cigar_op(cigar_element));
  }

  /**
   * @brief gets the length of an individual cigar element
   */
  inline static uint32_t cigar_oplen(const uint32_t cigar_element) {
    return bam_cigar_oplen(cigar_element);
  }

 private:
  std::shared_ptr<bam1_t> m_sam_record;   ///< sam record containing our cigar, potentially co-owned by multiple other objects
  uint32_t* m_cigar;                      ///< pointer to the start of the cigar in m_sam_record, cached for efficiency
  uint32_t m_num_cigar_elements;          ///< number of elements in our cigar

  static const char cigar_ops_as_chars[]; ///< static lookup table to convert CigarOperator enum values to chars.
};

}

#endif /* __gamgee__cigar__ */
