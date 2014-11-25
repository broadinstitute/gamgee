#ifndef gamgee__variant_filters_iterator__guard
#define gamgee__variant_filters_iterator__guard

#include "variant_filters.h"

#include "../utils/hts_memory.h"

#include "htslib/vcf.h"

#include <string>
#include <memory>

namespace gamgee {

/**
 * @brief simple random-access iterator class for VariantFilters objects
 */
class VariantFiltersIterator : public std::iterator<std::random_access_iterator_tag, std::string> {
 public:

  /**
   * @brief simple constructor used by the VariantFilters begin/end member functions
   * @param header an htslib variant header to keep shared ownership of the memory
   * @param body an htslib variant body to keep shared ownership of the memory
   * @param position current position in the iterator (starts at 0, normally)
   */
  VariantFiltersIterator(const std::shared_ptr<bcf_hdr_t>& header, const std::shared_ptr<bcf1_t>& body, const uint32_t position) : 
    m_header {header},
    m_body {body},
    m_position {position} 
  {}

  std::string operator[](int index) const { return utils::htslib_filter_name(m_header.get(), m_body.get(), index); }                                      ///< @brief returns a new string for the element at position 'index'
  std::string operator*() const { return utils::htslib_filter_name(m_header.get(), m_body.get(), m_position); }                                           ///< @brief returns a new string for the element at the current position
  VariantFiltersIterator operator++() { m_position++; return *this; }                                                                                     ///< @brief advances the iterator by one position
  VariantFiltersIterator operator--() { m_position--; return *this; }                                                                                     ///< @brief rewinds the iterator by one position
  VariantFiltersIterator operator+=(const int n) { m_position += n; return *this; }                                                                       ///< @brief advances the iterator by n positions
  VariantFiltersIterator operator-=(const int n) { m_position -= n; return *this; }                                                                       ///< @brief rewinds the iterator by n positions
  bool operator!=(const VariantFiltersIterator& other) { return m_position != other.m_position || m_body != other.m_body || m_header != other.m_header; } ///< @brief can be compared for equivalence using the equality/inequality operators
  bool operator==(const VariantFiltersIterator& other) { return m_position == other.m_position && m_body == other.m_body && m_header == other.m_header; } ///< @brief can be compared for equivalence using the equality/inequality operators
  bool operator<=(const VariantFiltersIterator& other) { return m_position <= other.m_position && m_body == other.m_body && m_header == other.m_header; } ///< @brief Can be compared with inequality relational operators (<, >, <= and >=).
  bool operator>=(const VariantFiltersIterator& other) { return m_position >= other.m_position && m_body == other.m_body && m_header == other.m_header; } ///< @brief Can be compared with inequality relational operators (<, >, <= and >=).
  bool operator< (const VariantFiltersIterator& other) { return m_position < other.m_position && m_body == other.m_body && m_header == other.m_header;  } ///< @brief Can be compared with inequality relational operators (<, >, <= and >=).
  bool operator> (const VariantFiltersIterator& other) { return m_position > other.m_position && m_body == other.m_body && m_header == other.m_header;  } ///< @brief Can be compared with inequality relational operators (<, >, <= and >=).

  uint32_t size() const { return uint32_t(m_body->d.n_flt); } ///< @brief returns the number of filters in this object

 private:
  std::shared_ptr<bcf_hdr_t> m_header; ///< shared ownership of the VariantHeader record memory so it stays alive while this object is in scope
  std::shared_ptr<bcf1_t> m_body;      ///< shared ownership of the Variant record memory so it stays alive while this object is in scope
  uint32_t m_position;                 ///< current position in the iterator
};

}

#endif // gamgee__variant_filters_iterator__guard
