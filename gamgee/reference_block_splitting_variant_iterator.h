#ifndef __gamgee__reference_block_splitting_variant_iterator__
#define __gamgee__reference_block_splitting_variant_iterator__

#include "variant.h"
#include "multiple_variant_iterator.h"

#include <list>

namespace gamgee {

/**
 * @brief Utility class to handle reference blocks while iterating over multiple variant files
 *
 * @warn This class is experimental/WIP
 */
class ReferenceBlockSplittingVariantIterator : public MultipleVariantIterator {
 public:

  /**
   * @brief creates an empty iterator (used for the end() method) 
   */
  ReferenceBlockSplittingVariantIterator() = default;

  /**
   * @brief initializes a new iterator based on a vector of input files (vcf or bcf)
   *
   * @param variant_files   vector of vcf/bcf files opened via the bcf_open() macro from htslib
   * @param variant_headers vector of variant headers corresponding to these files
   */
  ReferenceBlockSplittingVariantIterator(const std::vector<std::shared_ptr<htsFile>> variant_files, const std::vector<std::shared_ptr<bcf_hdr_t>> variant_headers);

  /**
   * @brief a ReferenceBlockSplittingVariantIterator move constructor guarantees all objects will have the same state.
   */
  ReferenceBlockSplittingVariantIterator(ReferenceBlockSplittingVariantIterator&&) = default;
  ReferenceBlockSplittingVariantIterator& operator=(ReferenceBlockSplittingVariantIterator&&) = default;

  /**
   * @brief a ReferenceBlockSplittingVariantIterator cannot be copied.
   */
  ReferenceBlockSplittingVariantIterator(const ReferenceBlockSplittingVariantIterator&) = delete;
  ReferenceBlockSplittingVariantIterator& operator=(const ReferenceBlockSplittingVariantIterator&) = delete;

  /**
   * @brief pseudo-inequality operator (needed by for-each loop)
   *
   * @warning this method does the minimal work necessary to determine that we have reached the end of iteration.
   * it is NOT a valid general-purpose inequality method.
   *
   * @param rhs the other ReferenceBlockSplittingVariantIterator to compare to
   *
   * @return whether both iterators have entered their end states
   */
  bool operator!=(const ReferenceBlockSplittingVariantIterator& rhs);

  /**
   * @brief dereference operator (needed by for-each loop)
   *
   * @return a reference to the iterator's Variant vector
   */
  std::vector<Variant>& operator*();

  /**
   * @brief advances the iterator, fetching the next vector
   *
   * @return a reference to the iterator's Variant vector
   */
  std::vector<Variant>& operator++();

 private:
  // fetches the next reference-block-split Variant vector
  // calls populate_pending() and populate_split_variants() as needed
  void fetch_next_split_vector();

  // populates the list of pending variants from the incoming vector of pre-split reference-block variants
  inline void populate_pending();

  // populates the vector of split variants from the list of pending variants, modifying the pending list as well
  inline void populate_split_variants();

  // holds the incoming reference-block variants before and during split operations
  std::vector<Variant> m_pending_variants;

  // caches next reference-block-split Variant vector
  std::vector<Variant> m_split_variants;

  unsigned int m_pending_chrom = UINT_MAX;
  unsigned int m_pending_start = UINT_MAX;
  unsigned int m_pending_min_end = UINT_MAX;
};

}  // end namespace gamgee

#endif	// __gamgee__reference_block_splitting_variant_iterator__
