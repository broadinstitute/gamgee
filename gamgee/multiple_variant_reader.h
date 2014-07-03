#ifndef gamgee__multiple_variant_reader__
#define gamgee__multiple_variant_reader__

#include "htslib/vcf.h"

#include "variant_header.h"
#include "utils/hts_memory.h"

namespace gamgee {

/**
 * @brief Utility class to read multiple VCF/BCF files with an appropriate iterator in a for-each loop.
 *
 * This class is designed to parse the files in for-each loops with the following signature:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (auto& record : MultipleVariantReader<MultipleVariantIterator>{filenames})
 *   do_something_with_record(record);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * You can also use it with the stdin or any other stream by using the default constructor 
 * or passing in an empty string for a filename, like so:
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (auto& record : MultipleVariantReader<MultipleVariantIterator>{filename1, stream2})
 *   do_something_with_record(record);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
template<class ITERATOR>
class MultipleVariantReader {
 public:

  /**
   * @brief a MultipleVariantReader cannot be copied safely, as it is iterating over streams.
   * @brief enables reading records in multiple files (vcf or bcf)
   *
   * @param filenames the names of the variant files
   */
  MultipleVariantReader(const std::vector<std::string>& filenames) :
    m_variant_files {},
    m_variant_headers {} {
      for (const auto& filename : filenames) {
    	auto* file_ptr = bcf_open(filename.empty() ? "-" : filename.c_str(), "r");
    	m_variant_files.push_back(file_ptr);
    	m_variant_headers.push_back(utils::make_shared_variant_header(bcf_hdr_read(file_ptr)));
      }
    }

  /**
   * @brief VariantReader should never be copied, but it can be moved around
   */
  MultipleVariantReader(MultipleVariantReader&& other) :
    m_variant_files {std::move(other.m_variant_files)},
    m_variant_headers {std::move(other.m_variant_headers)}
  {}

  /**
   * @brief closes the file streams if they are files
   */
  ~MultipleVariantReader() {
    for (const auto& file_ptr : m_variant_files)
      bcf_close(file_ptr);
  }

  /**
   * @brief a MultipleVariantReader cannot be copied safely, as it is iterating over streams.
   */
  MultipleVariantReader(const MultipleVariantReader&) = delete;

  /**
   * @brief creates an ITERATOR pointing at the start of the input streams (needed by for-each
   * loop)
   *
   * @return an ITERATOR ready to start parsing the files
   */
  ITERATOR begin() {
    return ITERATOR{m_variant_files, m_variant_headers};
  }

  /**
   * @brief creates an ITERATOR with nullified input streams (needed by for-each loop)
   *
   * @return an ITERATOR that will match the end status of the iterator at the end of the streams
   */
  ITERATOR end() {
    return ITERATOR{};
  }


  /**
   * @brief returns a variant header created from merging the headers of the files being read
   */
  // TODO temp passthrough
  inline VariantHeader header() { return VariantHeader{m_variant_headers[0]}; }

 private:
  std::vector<vcfFile*> m_variant_files;                            ///< vector of the internal file structures of the variant files
  std::vector<const std::shared_ptr<bcf_hdr_t> > m_variant_headers; ///< vector of the internal header structures of the variant files
};

}  // end namespace gamgee

#endif /* defined(__gamgee__multiple_variant_reader__) */
