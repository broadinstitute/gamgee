#ifndef __foghorn__sam_iterator__
#define __foghorn__sam_iterator__

#include <fstream>
#include "sam.h"
#include "htslib/sam.h"

namespace foghorn {

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
    SamIterator(samFile * sam_file_ptr, bam_hdr_t * sam_header_ptr);

    /**
     * @brief a SamIterator move constructor guarantees all objects will have the same state.
     */
    SamIterator(SamIterator&&);
    
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
     * @return a persistent Sam object independent from the iterator (a copy of the iterator's object)
     */
    Sam operator*();

    /**
     * @brief pre-fetches the next record and tests for end of file
     *
     * @return a reference to the object (it can be const& because this return value should only be used 
     *         by the for-each loop to check for the eof)
     */
    Sam operator++();

    /**
     * @brief takes care of all the memory allocations of the htslib sam reader interface
     */
    ~SamIterator();
    
  private:
    samFile * sam_file_ptr;     ///< pointer to the sam file
    bam_hdr_t * sam_header_ptr; ///< pointer to the sam header
    bam1_t * sam_record_ptr;    ///< pointer to the internal structure of the sam record. Useful to only allocate it once.
    Sam sam_record;             ///< temporary record to hold between fetch (operator++) and serve (operator*)
    
    Sam fetch_next_record();    ///< makes a new (through copy) Sam object that the user is free to use/keep without having to worry about memory management
};

}  // end namespace foghorn

#endif
