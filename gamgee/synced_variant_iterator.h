#ifndef gamgee__synced_variant_iterator__guard
#define gamgee__synced_variant_iterator__guard

#include "htslib/vcf.h"
#include "htslib/synced_bcf_reader.h"

#include "variant.h"

#include <memory>

namespace gamgee {

/**
 * @brief Utility class to enable for-each style iteration in the SyncedVariantReader class
 */
class SyncedVariantIterator {
 public:

  /**
   * @brief creates an empty iterator (used for the end() method) 
   */
  SyncedVariantIterator();

  /**
   * @brief initializes a new iterator based on a structure of synced vcf/bcf file readers
   *
   * @param synced_readers an htslib structure of synced vcf/bcf file readers created using the bcf_sr_init() function
    */
  SyncedVariantIterator(const std::shared_ptr<bcf_srs_t>& synced_readers);

  /**
   * @brief a SyncedVariantIterator should never be copied, but it can be moved around
   */
  SyncedVariantIterator(SyncedVariantIterator&& other) = default;
  SyncedVariantIterator& operator=(SyncedVariantIterator&& other) = default;

  /**
   * @brief a SyncedVariantIterator cannot be copied safely, as it is iterating over streams.
   */
  SyncedVariantIterator(const SyncedVariantIterator&) = delete;
  SyncedVariantIterator& operator=(const SyncedVariantIterator&) = delete;

  /**
   * @brief pseudo-inequality operator (needed by for-each loop)
   *
   * @warning this method does the minimal work necessary to determine that we have reached the end of iteration.
   * it is NOT a valid general-purpose inequality method.
   *
   * @param rhs the other SyncedVariantIterator to compare to
   *
   * @return whether both iterators have entered their end states
   */
  bool operator!=(const SyncedVariantIterator& rhs);

  /**
   * @brief dereference operator (needed by for-each loop)
   *
   * @return a persistent Variant vector independent from the iterator (a copy of the iterator's vector)
   */
  std::vector<Variant>& operator*();

  /**
   * @brief pre-fetches the next vector and tests for end of file
   *
   * @return a reference to the vector
   */
  std::vector<Variant>& operator++();

 private:
  std::shared_ptr<bcf_srs_t> m_synced_readers;                  ///< pointer to the synced readers of the variant files
  std::vector<Variant> m_variant_vector;                        ///< caches next Variant vector
  std::vector<std::shared_ptr<bcf_hdr_t>> m_headers_vector;     ///< caches each reader's htslib header

  void init_headers_vector();                                   ///< initializes m_variant_headers
  void fetch_next_record();                                     ///< fetches next Variant vector
};

}  // end namespace gamgee

#endif /* gamgee__synced_variant_iterator__guard */
