#ifndef gamgee__multiple_variant_reader__guard
#define gamgee__multiple_variant_reader__guard

#include "htslib/vcf.h"

#include "variant_header.h"
#include "utils/hts_memory.h"
#include "utils/variant_utils.h"

namespace gamgee {

/**
 * @brief Utility class to read multiple VCF/BCF files with an appropriate iterator in a for-each loop.
 *
 * @note result vectors will be ordered arbitrarily per location.  There is no guarantee or expectation that the vector
 * will be in the same order (or have the same number) as the inputs.  In the common case where all files have exactly
 * the same sites, it is recommended to use a zip iterator over SingleVariantReaders.  This will produce a deterministic
 * ordering and will be faster than MultipleVariantReader.
 *
 * This class is designed to parse files in for-each loops with the following signature:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (auto& vector : MultipleVariantReader<MultipleVariantIterator>{filenames})
 *   do_something_with_vector(vector);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * You can also use it with stdin or any other stream by using the default constructor
 * or passing in an empty string for a filename, like so:
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (auto& vector : MultipleVariantReader<MultipleVariantIterator>{filename1, stream2})
 *   do_something_with_vector(vector);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
template<class ITERATOR>
class MultipleVariantReader {
 public:

  /**
   * @brief enables reading records in multiple files (vcf or bcf)
   *
   * @param filenames the names of the variant files
   * @param validate_headers should we validate that the header files have identical chromosomes?  default = true
   */
  explicit MultipleVariantReader(const std::vector<std::string>& filenames, const bool validate_headers = true) :
    m_variant_files { },
    m_variant_header { }
  {
    init_reader(filenames, validate_headers);
  }

  /**
   * @brief enables reading records in multiple files (vcf or bcf)
   *
   * @param filenames the names of the variant files
   * @param validate_headers should we validate that the header files have identical chromosomes?  (must specify if using this constructor)
   * @param samples the list of samples you want included/excluded from your iteration
   * @param include whether you want these samples to be included or excluded from your iteration.  default = true (include)
   */
  MultipleVariantReader(const std::vector<std::string>& filenames, const bool validate_headers,
                        const std::vector<std::string>& samples, const bool include = true) :
    m_variant_files { },
    m_variant_header { nullptr }
  {
    init_reader(filenames, validate_headers);
    subset_variant_samples(m_variant_header.get(), samples, include);
  }

  /**
   * @brief helper function for constructors
   *
   * @param filenames the names of the variant files
   * @param validate_headers should we validate that the header files have identical chromosomes?
   */
  void init_reader(const std::vector<std::string>& filenames, const bool validate_headers) {
    for (const auto& filename : filenames) {
      // TODO? check for maximum one stream
      auto* file_ptr = bcf_open(filename.empty() ? "-" : filename.c_str(), "r");
      m_variant_files.push_back(utils::make_shared_hts_file(file_ptr));

      const auto& header_ptr = utils::make_shared_variant_header(bcf_hdr_read(file_ptr));
      if (!m_variant_header)
        m_variant_header = header_ptr;

      if (validate_headers)
        validate_header(header_ptr);
    }
  }

  /**
   * @brief MultipleVariantReader should never be copied, but it can be moved
   */
  MultipleVariantReader(MultipleVariantReader&& other) = default;
  MultipleVariantReader& operator=(MultipleVariantReader&& other) = default;

  /**
   * @brief a MultipleVariantReader cannot be copied safely, as it is iterating over streams.
   */
  MultipleVariantReader(const MultipleVariantReader&) = delete;
  MultipleVariantReader& operator=(const MultipleVariantReader& other) = delete;

  /**
   * @brief creates an ITERATOR pointing at the start of the input streams (needed by for-each
   * loop)
   *
   * @return an ITERATOR ready to start parsing the files
   */
  ITERATOR begin() const {
    return ITERATOR{m_variant_files, m_variant_header};
  }

  /**
   * @brief creates a default ITERATOR (needed by for-each loop)
   *
   * @return an ITERATOR that will match the end status of the iterator at the end of the streams
   */
  ITERATOR end() const {
    return ITERATOR{};
  }

  /**
   * @brief returns a representative header for the files being read
   *
   * @return a VariantHeader object constructed from the header of the first file read
   */
  const inline VariantHeader header() { return VariantHeader {m_variant_header}; }

  class HeaderException : public std::runtime_error {
   public:
    HeaderException() : std::runtime_error("Error: chromosomes in header files are inconsistent") { }
  };
 private:
  std::vector<std::shared_ptr<htsFile>> m_variant_files;        ///< vector of the internal file structures of the variant files
  std::shared_ptr<bcf_hdr_t> m_variant_header;                  ///< the internal header structure of the first variant file

  ///< confirms that the chromosomes in the headers of all of the input files are identical
  // TODO? only handles chromosome names, not lengths
  void validate_header(const std::shared_ptr<bcf_hdr_t>& other_header_ptr) {
    const auto& other_header = VariantHeader{other_header_ptr};
    if (header().chromosomes() != other_header.chromosomes())
      throw HeaderException{};
  }
};

}  // end namespace gamgee

#endif /* defined(gamgee__multiple_variant_reader__guard) */
