#ifndef gamgee__indexed_sam_reader__guard
#define gamgee__indexed_sam_reader__guard

#include "indexed_sam_iterator.h"
#include "utils/hts_memory.h"

#include "htslib/sam.h"

#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>

namespace gamgee {

/**
 * @brief Utility class to read a BAM/CRAM file with an appropriate Sam iterator from an indexed file
 * in a for-each loop. Intervals are passed in using a vector of string coordinates compatible with
 * Samtools. When iteration begins, the iterations (re-)starts at the beginning of the first interval.
 *
 * This class is designed to parse the file in for-each loops with the following signature:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (const auto& record : IndexedSamReader<IndexedSamIterator>{filename, interval_list})
 *   do_something_with_sam(record);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Most iterators have aliases defined by this module so you can use it like so:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (const auto& record : IndexedSingleSamReader{filename, interval_list})
 *   do_something_with_sam(record);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
template<class ITERATOR>
class IndexedSamReader {
  public:

    /**
     * @brief reads through all records in a file parsing them into Sam objects
     *
     * @param filename the name of the bam/cram file
     * @param interval_list Samtools style intervals to look for records
     */
    IndexedSamReader(const std::string& filename, std::vector<std::string> interval_list) :
      m_sam_file_ptr {sam_open(filename.c_str(), "r")}, // TODO: Not checking for errors.
      m_sam_index_ptr {sam_index_load(m_sam_file_ptr, filename.c_str())}, // TODO: Not checking for errors.
      m_sam_header_ptr { utils::make_shared_sam_header(sam_hdr_read(m_sam_file_ptr)) }, // TODO: Not checking for errors.
      m_interval_list {std::move(interval_list)} {
    }

    /**
     * @brief simple move assignment/construction
     *
     * we need a custom move constructor here because we need to set the file pointers to nullptr
     * otherwise the destructor will try to close the dangling pointers.
     *
     * @param other the other IndexedSamIterator to copy from
     */
    IndexedSamReader(IndexedSamReader&& other) :
      m_sam_file_ptr {std::move(other.m_sam_file_ptr)},
      m_sam_index_ptr {std::move(other.m_sam_index_ptr)},
      m_sam_header_ptr {std::move(other.m_sam_header_ptr)},
      m_interval_list {std::move(other.m_interval_list)} {
      other.m_sam_index_ptr = nullptr;
      other.m_sam_file_ptr = nullptr;
    }

    /**
     * @copydoc IndexedSamReader(IndexedSamReader&&)
     */
    IndexedSamReader& operator=(IndexedSamReader&& other) {
      m_sam_file_ptr = std::move(other.m_sam_file_ptr);
      m_sam_index_ptr = std::move(other.m_sam_index_ptr);
      m_sam_header_ptr = std::move(other.m_sam_header_ptr);
      m_interval_list = std::move(other.m_interval_list);
      other.m_sam_index_ptr = nullptr;
      other.m_sam_file_ptr = nullptr;
      return *this;
    }

    /**
     * @brief no copy construction/assignment allowed for iterators and readers
     */
    IndexedSamReader(IndexedSamReader& other) = delete;

    /**
     * @copydoc IndexedSamReader(IndexedSamReader&)
     */
    IndexedSamReader& operator=(IndexedSamReader& other) = delete;

    /**
     * @brief closes the file streams if there are some
     */
    ~IndexedSamReader() {
      if (m_sam_index_ptr)
        hts_idx_destroy(m_sam_index_ptr); // TODO: Copied from samtools. Is there a {bam,sam}_idx_destroy alias? Didn't find under "BAM/CRAM indexing" in sam.h
      if (m_sam_file_ptr)
        sam_close(m_sam_file_ptr);
    }

    /**
     * @brief creates a ITERATOR pointing at the start of the input stream (needed by for-each
     * loop)
     *
     * @return a ITERATOR ready to start parsing the file
     */
    ITERATOR begin() {
      if (m_interval_list.empty())
        return ITERATOR{};
      else
        return ITERATOR{m_sam_file_ptr, m_sam_index_ptr, m_sam_header_ptr, m_interval_list};
    }

    /**
     * @brief creates a ITERATOR with a nullified input stream (needed by for-each loop)
     *
     * @return a ITERATOR that will match the end status of the iterator at the end of the stream
     */
    ITERATOR end() {
      return ITERATOR{};
    }

    /**
     * @brief returns the header
     *
     * @return the header
     */
    inline SamHeader header() { return SamHeader{m_sam_header_ptr}; }

  private:
    samFile * m_sam_file_ptr;                    ///< pointer to the bam file
    hts_idx_t * m_sam_index_ptr;                 ///< pointer to the bam index
    std::shared_ptr<bam_hdr_t> m_sam_header_ptr; ///< pointer to the bam header
    std::vector<std::string> m_interval_list;    ///< intervals to iterate
};

using IndexedSingleSamReader = IndexedSamReader<IndexedSamIterator>;

}

#endif
