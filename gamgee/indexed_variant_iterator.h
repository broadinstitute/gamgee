#ifndef gamgee__indexed_variant_iterator__guard
#define gamgee__indexed_variant_iterator__guard

#include "variant_iterator.h"

#include "utils/hts_memory.h"

#include "htslib/vcf.h"

#include <memory>
#include <string>
#include <vector>

namespace gamgee {

class IndexedVariantIterator : public VariantIterator {
 public:

  static const std::vector<std::string> all_intervals;

  /**
   * @brief creates an empty iterator (used for the end() method)
   */
  IndexedVariantIterator();

  /**
   * @brief initializes a new iterator based on a file, an index, a header, and a vector of intervals
   *
   * @param file_ptr            shared pointer to a BCF file opened via the bcf_open() macro from htslib
   * @param index_ptr           shared pointer to a BCF file index (CSI) created with the bcf_index_load() macro from htslib
   * @param header_ptr          shared pointer to a BCF file header created with the bcf_hdr_read() macro from htslib
   * @param interval_list       vector of intervals represented by strings
   */
  IndexedVariantIterator(const std::shared_ptr<htsFile>& file_ptr,
                         const std::shared_ptr<hts_idx_t>& index_ptr,
                         const std::shared_ptr<bcf_hdr_t>& header_ptr,
                         const std::vector<std::string> interval_list = all_intervals);

  /**
   * @brief an IndexedVariantIterator cannot be copied safely, as it is iterating over a stream.
   */

  IndexedVariantIterator(const IndexedVariantIterator& other) = delete;
  IndexedVariantIterator& operator=(const IndexedVariantIterator& other) = delete;

  /**
   * @brief an IndexedVariantIterator can be moved
   */

  IndexedVariantIterator(IndexedVariantIterator&& other) = default;
  IndexedVariantIterator& operator=(IndexedVariantIterator&& other) = default;

  /**
   * @brief inequality operator (needed by for-each loop)
   *
   * @param rhs the other IndexedVariantIterator to compare to
   *
   * @return whether or not the two iterators are the same (e.g. have the same input file on the same
   * status and the same intervals)
   */
  bool operator!=(const IndexedVariantIterator& rhs);

 protected:
  void fetch_next_record() override;                                       ///< fetches next Variant record into existing htslib memory without making a copy

 private:
  std::shared_ptr<hts_idx_t> m_variant_index_ptr;                          ///< pointer to the internal structure of the index file
  std::vector<std::string> m_interval_list;                                ///< vector of intervals represented by strings
  std::vector<std::string>::const_iterator m_interval_iter;                ///< iterator for the interval list
  std::unique_ptr<hts_itr_t, utils::HtsIteratorDeleter> m_index_iter_ptr;  ///< pointer to the htslib BCF index iterator
};

}

#endif  /* defined(gamgee__indexed_variant_iterator__guard) */
