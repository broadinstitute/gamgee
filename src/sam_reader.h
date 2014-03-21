#ifndef foghorn__sam_reader__
#define foghorn__sam_reader__

#include <string>
#include <iostream>
#include <fstream>

#include "sam_iterator.h"
#include "htslib/sam.h"  // fix this ! 

namespace foghorn {

/**
 * @brief Utility class to read a SAM/BAM/CRAM file with an appropriate Sam iterator from a stream 
 * (e.g. file, stdin, ...) in a for-each loop.
 *
 * This class is designed to parse the file in for-each loops with the following signature:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (auto record : SamReader<SamIterator>(filename))
 *   do_something_with_sam(record);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * You can also use it with the stdin or any other stream by using the default constructor 
 * or passing in an empty string for a filename, like so:
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (auto pair : SamReader<SamPairIterator>())
 *   do_something_with_pair(pair);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
template<class ITERATOR>
class SamReader {
  public:

    /**
     * @brief reads through all records in a file ( or sam) parsing them into Sam
     * objects
     *
     * @param filename the name of the sam file
     */
    SamReader(const std::string& filename) :
      sam_file_ptr {sam_open(filename.empty() ? "-" : filename.c_str(), "r")},
      sam_header_ptr {sam_hdr_read(sam_file_ptr)}
    {}

    /**
     * @brief closes the file stream if there is one (in case we are reading a sam file)
     */
    ~SamReader() {
      sam_close(sam_file_ptr);
      bam_hdr_destroy(sam_header_ptr);
    }

    /**
     * @brief a SamReader cannot be copied safely, as it is iterating over a stream.
     */
    SamReader(const SamReader&) = delete;

    /**
     * @brief creates a ITERATOR pointing at the start of the input stream (needed by for-each
     * loop)
     *
     * @return a ITERATOR ready to start parsing the file
     */
    ITERATOR begin() {
      return ITERATOR{sam_file_ptr, sam_header_ptr};
    }

    /**
     * @brief creates a ITERATOR with a nullified input stream (needed by for-each loop)
     *
     * @return a ITERATOR that will match the end status of the iterator at the end of the stream
     */
    ITERATOR end() {
      return ITERATOR{};
    }

  private:
    samFile *sam_file_ptr;     ///< pointer to the internal file structure of the sam/bam/cram file
    bam_hdr_t *sam_header_ptr; ///< pointer to the internal header structure of the sam/bam/cram file
};

}  // end of namespace

#endif /* defined(__foghorn__sam_iterator__) */
