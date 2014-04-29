#ifndef __gamgee__cigar__
#define __gamgee__cigar__

#include "htslib/sam.h"

#include <memory>
#include <string>

namespace gamgee {

// Order of the operators in this enum must match the order in BAM_CIGAR_STR from htslib/sam.h
enum class CigarOperator { M, I, D, N, S, H, P, EQ, X, B };

class Cigar {
public:
  explicit Cigar(const std::shared_ptr<bam1_t>& sam_record);
  Cigar(const Cigar& other);
  Cigar(Cigar&& other);
  Cigar& operator=(const Cigar& other);
  Cigar& operator=(Cigar&& other);

  // Default destruction is sufficient, since our shared_ptr will handle deallocation
  ~Cigar() = default;

  uint32_t operator[](const uint32_t index) const;
  uint32_t size() const { return m_num_cigar_elements; }
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
  // Sam record containing our cigar, potentially co-owned by multiple other objects
  std::shared_ptr<bam1_t> m_sam_record;

  // Pointer to the start of the cigar in m_sam_record, cached for efficiency
  uint32_t* m_cigar;

  // Number of elements in our cigar
  uint32_t m_num_cigar_elements;

  static const char cigar_ops_as_chars[];
};

}

#endif /* __gamgee__cigar__ */
