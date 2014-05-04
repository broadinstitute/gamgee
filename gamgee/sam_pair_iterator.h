#ifndef __gamgee__sam_pair_iterator__
#define __gamgee__sam_pair_iterator__

#include "sam.h"

#include "htslib/sam.h"

#include <fstream>
#include <queue>
#include <memory>

namespace gamgee {

/**
 * @brief Utility class to enable for-each style iteration by pairs in the SamReader class
 */
class SamPairIterator {
  public:

    /**
     * @brief creates an empty iterator (used for the end() method) 
     */
    SamPairIterator();

    /**
     * @brief initializes a new iterator based on an input stream (e.g. sam/a file, stdin, ...)
     *
     * @param sam_file_ptr   pointer to a sam file opened via the sam_open() macro from htslib
     * @param sam_header_ptr pointer to a sam file header created with the sam_hdr_read() macro from htslib
     */
    SamPairIterator(samFile * sam_file_ptr, const std::shared_ptr<bam_hdr_t>& sam_header_ptr);

    /**
     * @brief a SamPairIterator move constructor guarantees all objects will have the same state.
     */
    SamPairIterator(SamPairIterator&& other);
    
    /**
     * @brief inequality operator (needed by for-each loop)
     *
     * @param rhs the other SamPairIterator to compare to
     *
     * @return whether or not the two iterators are the same (e.g. have the same input stream on the same
     * status)
     */
    bool operator!=(const SamPairIterator& rhs);

    /**
     * @brief dereference operator (needed by for-each loop)
     *
     * @return a persistent Sam object independent from the iterator (a copy of the iterator's object)
     */
    std::pair<Sam,Sam> operator*();

    /**
     * @brief pre-fetches the next record and tests for end of file
     *
     * @return a reference to the object (it can be const& because this return value should only be used 
     *         by the for-each loop to check for the eof)
     */
    std::pair<Sam,Sam> operator++();

  private:
    using SamPtrQueue = std::queue<std::shared_ptr<bam1_t>>;

    SamPtrQueue m_supp_alignments;                     ///< queue to hold the supplementary alignments temporarily while processing the pairs
    samFile * m_sam_file_ptr;                          ///< pointer to the sam file
    const std::shared_ptr<bam_hdr_t> m_sam_header_ptr; ///< pointer to the sam header
    std::shared_ptr<bam1_t> m_sam_record_ptr1;         ///< pointer to the internal structure of the sam record. Useful to only allocate it once.
    std::shared_ptr<bam1_t> m_sam_record_ptr2;         ///< pointer to the internal structure of the sam record. Useful to only allocate it once.
    std::pair<Sam,Sam> m_sam_records;                  ///< temporary record to hold between fetch (operator++) and serve (operator*)

    std::pair<Sam,Sam> fetch_next_pair();              ///< makes a new (through copy) pair of Sam objects that the user is free to use/keep without having to worry about memory management
    bool read_sam(std::shared_ptr<bam1_t>& record_ptr);                 ///< reads a sam record and checks for the end-of-file invalidating the file and header pointers if necessary
    Sam make_sam(std::shared_ptr<bam1_t>& record_ptr);                  ///< creates a sam record from the internal data
    Sam next_primary_alignment(std::shared_ptr<bam1_t>& record_ptr);
    std::pair<Sam,Sam> next_supplementary_alignment();
};

}  // end namespace gamgee

#endif
