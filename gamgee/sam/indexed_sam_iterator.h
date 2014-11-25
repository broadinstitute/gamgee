#ifndef gamgee__indexed_sam_iterator__guard
#define gamgee__indexed_sam_iterator__guard

#include "sam.h"

#include "../utils/hts_memory.h"

#include "htslib/sam.h"

#include <memory>
#include <string>
#include <vector>

namespace gamgee {

/**
 * @brief Utility class to enable for-each style iteration in the IndexedSamReader class
 */
class IndexedSamIterator {
  public:
    /**
     * @brief creates an empty iterator (used for the end() method)
     */
    IndexedSamIterator();

    /**
     * @brief initializes a new iterator based on an input file stream
     *
     * @param sam_file_ptr   pointer to a bam/cram file opened via the sam_open() macro from htslib
     * @param sam_index_ptr  pointer to a bam/cram file opened via the sam_index_load() macro from htslib
     * @param sam_header_ptr pointer to a bam/cram file header created with the sam_hdr_read() macro from htslib
     * @param interval_list  vector of intervals compatible with sam_itr_querys, hts_parse_reg, etc.
     */
    IndexedSamIterator(const std::shared_ptr<htsFile>& sam_file_ptr, const std::shared_ptr<hts_idx_t>& sam_index_ptr,
        const std::shared_ptr<bam_hdr_t>& sam_header_ptr, const std::vector<std::string>& interval_list);

    /**
     * @brief iterators and readers can be moved
     */
    IndexedSamIterator(IndexedSamIterator&& other) = default;
    IndexedSamIterator& operator=(IndexedSamIterator&& other) = default;

    /**
     * @brief no copy construction/assignment allowed for iterators and readers
     */
    IndexedSamIterator(const IndexedSamIterator& other) = delete;
    IndexedSamIterator& operator=(const IndexedSamIterator& other) = delete;

    /**
     * @brief inequality operator (needed by for-each loop)
     *
     * @param rhs the other SamIterator to compare to
     *
     * @return whether or not the two iterators are the same (e.g. have the same input stream on the same
     * status)
     */
    bool operator!=(const IndexedSamIterator& rhs);

    /**
     * @brief dereference operator (needed by for-each loop)
     *
     * @return a Sam object by reference, valid until the next record is fetched (the iterator re-uses memory at each iteration)
     */
    Sam& operator*();

    /**
     * @brief pre-fetches the next record and tests for end of file
     *
     * @return a reference to the object (it can be const& because this return value should only be used
     *         by the for-each loop to check for the eof)
     */
    Sam& operator++();

    const std::string& current_interval() const;

  private:
    std::shared_ptr<htsFile> m_sam_file_ptr;                ///< pointer to the bam file
    std::shared_ptr<hts_idx_t> m_sam_index_ptr;             ///< pointer to the bam index
    std::shared_ptr<bam_hdr_t> m_sam_header_ptr;            ///< pointer to the bam header
    std::vector<std::string> m_interval_list;               ///< intervals to iterate
    std::vector<std::string>::iterator m_interval_iterator; ///< temporary interval to hold between sam_itr_querys and serve fetch_next_record
    std::unique_ptr<hts_itr_t, utils::HtsIteratorDeleter> m_sam_itr_ptr; ///< temporary iterator to hold between sam_itr_querys and serve fetch_next_record
    std::shared_ptr<bam1_t> m_sam_record_ptr;               ///< pointer to the internal structure of the sam record. Useful to only allocate it once.
    Sam m_sam_record;                                       ///< temporary record to hold between fetch (operator++) and serve (operator*)

    void fetch_next_record();                               ///< fetches next Sam record into existing htslib memory without making a copy
};

}

#endif
