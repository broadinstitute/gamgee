#ifndef gamgee__indexed_variant_reader__guard
#define gamgee__indexed_variant_reader__guard

#include "indexed_variant_iterator.h"
#include "utils/hts_memory.h"

#include "htslib/vcf.h"

#include <memory>
#include <string>
#include <vector>

namespace gamgee {

  /**
   * @brief Utility class to read an indexed BCF file by intervals using an appropriate Variant iterator
   * in a for-each loop.
   *
   * NOTE: this will only parse BCF files with CSI indices
   *
   * This class is designed to parse the file in for-each loops with the following signature:
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * for (auto& record : IndexedVariantReader<IndexedVariantIterator>{filename, intervals})
   *   do_something_with_record(record);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
template<class ITERATOR>
class IndexedVariantReader {
 public:

  /**
   * @brief reads through all records in a file matching one of the given intervals,
   * parsing them into Variant objects
   *
   * @param filename the name of the variant file
   * @param interval_list a vector of intervals represented by strings.  Empty vector for all intervals.
   *
   */
  IndexedVariantReader(const std::string& filename, const std::vector<std::string>& interval_list) :
    m_variant_file_ptr { utils::make_shared_hts_file(bcf_open(filename.c_str(), "r")) },
    m_variant_index_ptr { utils::make_shared_hts_index(bcf_index_load(filename.c_str())) },
    m_variant_header_ptr { utils::make_shared_variant_header(bcf_hdr_read(m_variant_file_ptr.get())) },
    m_interval_list { interval_list }
  {}

  /**
   * @brief an IndexedVariantReader cannot be copied safely, as it is iterating over a stream.
   */

  IndexedVariantReader(const IndexedVariantReader& other) = delete;
  IndexedVariantReader& operator=(const IndexedVariantReader& other) = delete;

  /**
   * @brief an IndexedVariantReader can be moved
   */

  IndexedVariantReader(IndexedVariantReader&& other) = default;
  IndexedVariantReader& operator=(IndexedVariantReader&& other) = default;

  ITERATOR begin() const {
    return ITERATOR{ m_variant_file_ptr, m_variant_index_ptr, m_variant_header_ptr, m_interval_list };
  }

  ITERATOR end() const {
    return ITERATOR{};
  }

  /**
   * @brief returns the variant header of the file being read
   */
  inline VariantHeader header() const { return VariantHeader{m_variant_header_ptr}; }

 private:
  std::shared_ptr<vcfFile> m_variant_file_ptr;        ///< pointer to the internal structure of the variant file
  std::shared_ptr<hts_idx_t> m_variant_index_ptr;     ///< pointer to the internal structure of the index file
  std::shared_ptr<bcf_hdr_t> m_variant_header_ptr;    ///< pointer to the internal structure of the header file
  std::vector<std::string> m_interval_list;           ///< vector of intervals represented by strings
};

}

#endif  /* defined(gamgee__indexed_variant_reader__guard) */
