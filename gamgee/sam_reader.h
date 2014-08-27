#ifndef gamgee__sam_reader__guard
#define gamgee__sam_reader__guard

#include "sam_iterator.h"
#include "sam_pair_iterator.h"
#include "utils/hts_memory.h"
#include "exceptions.h"

#include "htslib/sam.h"

#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>

namespace gamgee {

/**
 * @brief Utility class to read a SAM/BAM/CRAM file with an appropriate Sam iterator from a stream 
 * (e.g. file, stdin, ...) in a for-each loop.
 *
 * This class is designed to parse the file in for-each loops with the following signature:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (auto& record : SamReader<SamIterator>{filename})
 *   do_something_with_sam(record);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * You can also use it with the stdin or any other stream by using the default constructor 
 * or passing in an empty string for a filename, like so:
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (auto& pair : SamReader<SamPairIterator>{filename})
 *   do_something_with_pair(pair);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Most iterators have aliases defined by this module so you can use it like so:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (auto& pair : SingleSamReader{filename})
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
      m_sam_file_ptr {},
      m_sam_header_ptr {}
    {
      init_reader(filename);
    }

    /**
     * @brief reads through all records in a file ( or sam) parsing them into Sam
     * objects
     *
     * @param filenames a vector containing a single element: the name of the sam file
     */
    SamReader(const std::vector<std::string>& filenames) :
      m_sam_file_ptr {},
      m_sam_header_ptr {}
    {
      if (filenames.size() > 1)
        throw SingleInputException{"filenames", filenames.size()};
      if (!filenames.empty())
        init_reader(filenames.front());
    }

    /**
     * @brief no copy construction/assignment allowed for iterators and readers
     */
    SamReader(const SamReader& other) = delete;
    SamReader& operator=(const SamIterator&) = delete;

    /**
     * @brief a SamReader move constructor guarantees all objects will have the same state.
     */
    SamReader(SamReader&&) = default;
    SamReader& operator=(SamReader&&) = default;

    /**
     * @brief creates a ITERATOR pointing at the start of the input stream (needed by for-each
     * loop)
     *
     * @return a ITERATOR ready to start parsing the file
     */
    ITERATOR begin() {
      return ITERATOR{m_sam_file_ptr, m_sam_header_ptr};
    }

    /**
     * @brief creates a ITERATOR with a nullified input stream (needed by for-each loop)
     *
     * @return a ITERATOR that will match the end status of the iterator at the end of the stream
     */
    ITERATOR end() {
      return ITERATOR{};
    }

    inline SamHeader header() { return SamHeader{m_sam_header_ptr}; }

  private:
    std::shared_ptr<htsFile> m_sam_file_ptr;     ///< pointer to the internal file structure of the sam/bam/cram file
    std::shared_ptr<bam_hdr_t> m_sam_header_ptr; ///< pointer to the internal header structure of the sam/bam/cram file

    /**
     * @brief initialize the SamReader (helper function for constructors)
     *
     * @param filename the name of the variant file
     */
    void init_reader (const std::string& filename) {
      auto* file_ptr = sam_open(filename.empty() ? "-" : filename.c_str(), "r");
      m_sam_file_ptr  = utils::make_shared_hts_file(file_ptr);
      m_sam_header_ptr = utils::make_shared_sam_header(sam_hdr_read(file_ptr));
    }
};

using SingleSamReader = SamReader<SamIterator>;
using PairSamReader = SamReader<SamPairIterator>;

}  // end of namespace

#endif /* defined(gamgee__sam_reader__guard) */
