#ifndef __gamgee__multiple_variant_iterator__
#define __gamgee__multiple_variant_iterator__

#include "htslib/vcf.h"

#include "variant.h"
#include "variant_iterator.h"

#include <memory>
#include <queue>

namespace gamgee {

/**
 * @brief Utility class to enable for-each style iteration in the MultipleVariantReader class
 */
class MultipleVariantIterator {
 public:

  /**
   * @brief creates an empty iterator (used for the end() method) 
   */
  MultipleVariantIterator();

  /**
   * @brief initializes a new iterator based on a vector if input files (vcf or bcf)
   *
   * @param variant_files   vector of vcf/bcf files opened via the bcf_open() macro from htslib
   * @param variant_headers vector of shared pointers to vcf/bcf file headers created with the bcf_hdr_read() macro from htslib
   */
  MultipleVariantIterator(const std::vector<vcfFile*> variant_files, const std::vector<const std::shared_ptr<bcf_hdr_t> > variant_headers);

  /**
   * @brief a MultipleVariantIterator move constructor guarantees all objects will have the same state.
   */
  MultipleVariantIterator(MultipleVariantIterator&&);

  /**
   * @brief inequality operator (needed by for-each loop)
   *
   * @param rhs the other MultipleVariantIterator to compare to
   *
   * @return whether or not the two iterators are the same (e.g. have the same input streams on the same
   * status)
   */
  bool operator!=(const MultipleVariantIterator& rhs);

  /**
   * @brief dereference operator (needed by for-each loop)
   *
   * @return a persistent Variant vector independent from the iterator (a copy of the iterator's vector)
   */
  std::vector<Variant>& operator*();

  /**
   * @brief pre-fetches the next vector and tests for end of file
   *
   * @return a reference to the vector (it can be const& because this return value should only be used
   *         by the for-each loop to check for the eof)
   */
  std::vector<Variant>& operator++();

 private:
  std::vector<Variant> fetch_next_record();				///< fetches next Variant vector
  std::priority_queue<std::shared_ptr<VariantIterator> > m_queue;	///< the individual file iterators
  std::vector<Variant> m_variant_vector;				///< caches next Variant vector
};

}  // end namespace gamgee

#endif	// __gamgee__multiple_variant_iterator__
