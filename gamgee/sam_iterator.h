#ifndef gamgee__sam_iterator__guard
#define gamgee__sam_iterator__guard

#include "sam.h"

#include "htslib/sam.h"

#include <memory>

namespace gamgee {

/**
 * @brief Utility class to enable for-each style iteration in the SamReader class
 */
class SamIterator {
  public:

    /**
     * @brief creates an empty iterator (used for the end() method) 
     */
    SamIterator();

    /**
     * @brief initializes a new iterator based on an input stream (e.g. sam/a file, stdin, ...)
     *
     * @param sam_file_ptr   pointer to a sam file opened via the sam_open() macro from htslib
     * @param sam_header_ptr pointer to a sam file header created with the sam_hdr_read() macro from htslib
     */
    SamIterator(const std::shared_ptr<htsFile>& sam_file_ptr, const std::shared_ptr<bam_hdr_t>& sam_header_ptr);

    /**
     * @brief no copy construction/assignment allowed for readers and iterators
     */
    SamIterator(const SamIterator&) = delete;
    SamIterator& operator=(const SamIterator&) = delete;

    /**
     * @brief a SamIterator move constructor guarantees all objects will have the same state.
     */
    SamIterator(SamIterator&&) = default;
    SamIterator& operator=(SamIterator&&) = default;
    
    /**
     * @brief inequality operator (needed by for-each loop)
     *
     * @param rhs the other SamIterator to compare to
     *
     * @return whether or not the two iterators are the same (e.g. have the same input stream on the same
     * status)
     */
    bool operator!=(const SamIterator& rhs);

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

  private:
    std::shared_ptr<htsFile> m_sam_file_ptr;     ///< pointer to the sam file
    std::shared_ptr<bam_hdr_t> m_sam_header_ptr; ///< pointer to the sam header
    std::shared_ptr<bam1_t> m_sam_record_ptr;    ///< pointer to the internal structure of the sam record. Useful to only allocate it once.
    Sam m_sam_record;                            ///< temporary record to hold between fetch (operator++) and serve (operator*)

    void fetch_next_record();                    ///< fetches next Sam record into existing htslib memory without making a copy
};

}  // end namespace gamgee

#endif // gamgee__sam_iterator__guard
