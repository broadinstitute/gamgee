#ifndef __gamgee__variant_iterator__
#define __gamgee__variant_iterator__

#include "variant.h"

#include "htslib/vcf.h"

#include <memory>

namespace gamgee {

/**
 * @brief Utility class to enable for-each style iteration in the VariantReader class
 */
class VariantIterator {
 public:

  /**
   * @brief creates an empty iterator (used for the end() method) 
   */
  VariantIterator();

  /**
   * @brief initializes a new iterator based on an input stream (e.g. a vcf/bcf file, stdin, ...)
   *
   * @param variant_file_ptr   pointer to a vcf/bcf file opened via the bcf_open() macro from htslib
   * @param variant_header_ptr shared pointer to a vcf/bcf file header created with the bcf_hdr_read() macro from htslib
   */
  VariantIterator(vcfFile* variant_file_ptr, const std::shared_ptr<bcf_hdr_t>& variant_header_ptr);

  /**
   * @brief a VariantIterator move constructor guarantees all objects will have the same state.
   */
  VariantIterator(VariantIterator&&);

  /**
   * @brief inequality operator (needed by for-each loop)
   *
   * @param rhs the other VariantIterator to compare to
   *
   * @return whether or not the two iterators are the same (e.g. have the same input stream on the same
   * status)
   */
  bool operator!=(const VariantIterator& rhs);

  /**
   * @brief dereference operator (needed by for-each loop)
   *
   * @return a persistent Variant object independent from the iterator (a copy of the iterator's object)
   */
  Variant& operator*();

  /**
   * @brief pre-fetches the next record and tests for end of file
   *
   * @return a reference to the object (it can be const& because this return value should only be used 
   *         by the for-each loop to check for the eof)
   */
  Variant& operator++();

 private:
  vcfFile * m_variant_file_ptr;                    ///< pointer to the vcf/bcf file
  std::shared_ptr<bcf_hdr_t> m_variant_header_ptr; ///< pointer to the variant header
  std::shared_ptr<bcf1_t> m_variant_record_ptr;    ///< pointer to the internal structure of the variant record. Useful to only allocate it once.
  Variant m_variant_record;                        ///< temporary record to hold between fetch (operator++) and serve (operator*)

  Variant fetch_next_record();                     ///< fetches next Variant record into existing htslib memory without making a copy
};

}  // end namespace gamgee

#endif
