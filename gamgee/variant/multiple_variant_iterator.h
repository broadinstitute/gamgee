#ifndef gamgee__multiple_variant_iterator__guard
#define gamgee__multiple_variant_iterator__guard

#include "htslib/vcf.h"

#include "variant.h"
#include "variant_iterator.h"

#include <memory>
#include <queue>

namespace gamgee {

using VariantIteratorIndexPair = std::pair<std::shared_ptr<VariantIterator>, uint32_t>;
using VariantIndexPair = std::pair<Variant, uint32_t>;

/**
 * @brief Utility class to enable for-each style iteration in the MultipleVariantReader class
 */
class MultipleVariantIterator {
 public:

  /**
   * @brief creates an empty iterator (used for the end() method) 
   */
  MultipleVariantIterator() = default;

  /**
   * @brief initializes a new iterator based on a vector of input files (vcf or bcf)
   *
   * @param variant_files   vector of vcf/bcf files opened via the bcf_open() macro from htslib
   * @param variant_headers vector of headers corresponding to the files
   */
  MultipleVariantIterator(const std::vector<std::shared_ptr<htsFile>>& variant_files, const std::vector<std::shared_ptr<bcf_hdr_t>>& variant_headers);

  /**
   * @brief a MultipleVariantIterator move constructor guarantees all objects will have the same state.
   */
  MultipleVariantIterator(MultipleVariantIterator&&) = default;
  MultipleVariantIterator& operator=(MultipleVariantIterator&& other) = default;

  /**
   * @brief a MultipleVariantIterator cannot be copied safely, as it is iterating over streams.
   */
  MultipleVariantIterator(const MultipleVariantIterator&) = delete;
  MultipleVariantIterator& operator=(const MultipleVariantIterator& other) = delete;

  /**
   * @brief pseudo-inequality operator (needed by for-each loop)
   *
   * @warning this method does the minimal work necessary to determine that we have reached the end of iteration.
   * it is NOT a valid general-purpose inequality method.
   *
   * @param rhs the other MultipleVariantIterator to compare to
   *
   * @return whether both iterators have entered their end states
   */
  bool operator!=(const MultipleVariantIterator& rhs);

  /**
   * @brief dereference operator (needed by for-each loop)
   *
   * @return a reference to the iterator's Variant vector
   */
  std::vector<VariantIndexPair>& operator*();

  /**
   * @brief advances the iterator, fetching the next vector
   *
   * @return a reference to the iterator's Variant vector
   */
  std::vector<VariantIndexPair>& operator++();

 private:
  // fetches the next Variant vector
  void fetch_next_vector();

  // comparison class for genomic locations in the priority queue
  class Comparator {
   public:
    bool operator()(const VariantIteratorIndexPair& left, const VariantIteratorIndexPair& right);
  };

  // the individual file iterators
  std::priority_queue<VariantIteratorIndexPair, std::vector<VariantIteratorIndexPair>, Comparator> m_queue;

  // caches next Variant vector
  std::vector<VariantIndexPair> m_variant_vector;
};

}  // end namespace gamgee

#endif	// gamgee__multiple_variant_iterator__guard
