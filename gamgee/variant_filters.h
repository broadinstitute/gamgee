#ifndef gamgee__variant_filters__guard
#define gamgee__variant_filters__guard

#include "variant_filters_iterator.h"
#include "utils/hts_memory.h"

#include "htslib/vcf.h"

#include <string>
#include <memory>

namespace gamgee {

/**
 * @brief class to manipulate filter field objects without making copies. 
 *
 * It provides all the functionality that a vector of strings should provide.
 * You can do random access or iteration style access on it. The iterator is
 * implemented as a random iterator for maximum flexibility with standard
 * library functions
 */
class VariantFilters {
 public:
  /**
   * @brief standard constructor used by the Variant API
   * @param header an htslib variant header to keep shared ownership of the memory
   * @param body an htslib variant body to keep shared ownership of the memory
   */
  VariantFilters(const std::shared_ptr<bcf_hdr_t>& header, const std::shared_ptr<bcf1_t>& body) : m_header {header}, m_body {body} {} 

  /**
   * @brief random access operator
   * @param index the filter index 
   */
  std::string operator[](int index) const { return utils::htslib_filter_name(m_header.get(), m_body.get(), index); }

  /**
   * @brief returns the number of filters in the filter field
   */
  uint32_t size() const { return uint32_t(m_body->d.n_flt); }

  /**
   * @brief returns an iterator pointing to the first element in the list of filters.
   */
  VariantFiltersIterator begin() const {return VariantFiltersIterator{m_header, m_body, 0}; }

  /**
   * @brief Returns an iterator referring to one-past-the-last element in the list of filters.
   */
  VariantFiltersIterator end() const {return VariantFiltersIterator{m_header, m_body, uint32_t(m_body->d.n_flt)};}

 private:
  std::shared_ptr<bcf_hdr_t> m_header; ///< shared ownership of the VariantHeader record memory so it stays alive while this object is in scope
  std::shared_ptr<bcf1_t> m_body;      ///< shared ownership of the Variant record memory so it stays alive while this object is in scope

};


}

#endif // gamgee__variant_filters__guard
