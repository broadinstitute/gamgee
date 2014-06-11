#ifndef gamgee__variant_reader__
#define gamgee__variant_reader__

#include "variant_header.h"
#include "variant_iterator.h"
#include "utils/hts_memory.h"

#include "htslib/vcf.h"

#include <string>
#include <iostream>
#include <fstream>
#include <memory>


namespace gamgee {

/**
 * @brief Utility class to read a VCF/BCF file with an appropriate Variant iterator from a stream 
 * (e.g. file, stdin, ...) in a for-each loop.
 *
 * This class is designed to parse the file in for-each loops with the following signature:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (auto& record : VariantReader<VariantIterator>{filename})
 *   do_something_with_record(record);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * You can also use it with the stdin or any other stream by using the default constructor 
 * or passing in an empty string for a filename, like so:
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (auto& record : VariantReader<VariantIterator>{filename})
 *   do_something_with_record(record);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Most iterators have aliases defined by this module so you can use it like so:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (auto& record : SingleVariantReader{filename})
 *   do_something_with_record(record);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
template<class ITERATOR>
class VariantReader {
 public:

  /**
   * @brief reads through all records in a file (vcf or bcf) parsing them into Variant
   * objects
   *
   * @param filename the name of the variant file
   */
  VariantReader(const std::string& filename) :
    m_variant_file_ptr {bcf_open(filename.empty() ? "-" : filename.c_str(), "r")},
    m_variant_header_ptr { utils::make_shared_variant_header(bcf_hdr_read(m_variant_file_ptr)) }
  {}

  /**
   * @brief VariantReader should never be copied, but it can be moved around
   */
  VariantReader(VariantReader&& other) :
    m_variant_file_ptr {std::move(other.m_variant_file_ptr)},
    m_variant_header_ptr {std::move(other.m_variant_header_ptr)}
  {}

  /**
   * @brief closes the file stream if there is one (in case we are reading a variant file)
   */
  ~VariantReader() {
    bcf_close(m_variant_file_ptr);
  }

  /**
   * @brief a VariantReader cannot be copied safely, as it is iterating over a stream.
   */
  VariantReader(const VariantReader&) = delete;

  /**
   * @brief creates a ITERATOR pointing at the start of the input stream (needed by for-each
   * loop)
   *
   * @return a ITERATOR ready to start parsing the file
   */
  ITERATOR begin() const {
    return ITERATOR{m_variant_file_ptr, m_variant_header_ptr};
  }

  /**
   * @brief creates a ITERATOR with a nullified input stream (needed by for-each loop)
   *
   * @return a ITERATOR that will match the end status of the iterator at the end of the stream
   */
  ITERATOR end() const {
    return ITERATOR{};
  }

  /**
   * @brief returns the variant header of the file being read
   */
  inline VariantHeader header() const { return VariantHeader{m_variant_header_ptr}; }

 private:
  vcfFile* m_variant_file_ptr;                           ///< pointer to the internal file structure of the variant/bam/cram file
  const std::shared_ptr<bcf_hdr_t> m_variant_header_ptr; ///< pointer to the internal header structure of the variant/bam/cram file
};

using SingleVariantReader = VariantReader<VariantIterator>;

}  // end of namespace

#endif /* defined(__gamgee__variant_reader__) */
