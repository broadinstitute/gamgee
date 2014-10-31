#ifndef gamgee__multiple_variant_reader__guard
#define gamgee__multiple_variant_reader__guard

#include "htslib/vcf.h"

#include "variant_header.h"
#include "utils/hts_memory.h"
#include "utils/variant_utils.h"
#include "exceptions.h"

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
    m_variant_headers { },
    m_combined_header { }
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
    m_variant_headers { },
    m_combined_header { }
  {
    init_reader(filenames, validate_headers);
    subset_variant_samples(m_combined_header.get(), samples, include);
  }

  /**
   * @brief helper function for constructors
   *
   * @param filenames the names of the variant files
   * @param validate_headers should we validate that the header files have identical chromosomes?
   */
  void init_reader(const std::vector<std::string>& filenames, const bool validate_headers) {
    m_variant_files.reserve(filenames.size());
    m_variant_headers.reserve(filenames.size());

    for (const auto& filename : filenames) {
      // TODO? check for maximum one stream
      auto* file_ptr = bcf_open(filename.empty() ? "-" : filename.c_str(), "r");
      if ( file_ptr == nullptr ) {
        throw FileOpenException{filename};
      }
      m_variant_files.push_back(utils::make_shared_hts_file(file_ptr));

      auto* header_raw_ptr = bcf_hdr_read(file_ptr);
      if ( header_raw_ptr == nullptr ) {
        throw HeaderReadException{filename};
      }
      const auto& header_ptr = utils::make_shared_variant_header(header_raw_ptr);
      m_variant_headers.push_back(header_ptr);

      if (m_combined_header) {
        if (validate_headers)
          validate_header(header_ptr);

        merge_variant_headers(m_combined_header, header_ptr);
      }
      else
        m_combined_header = utils::make_shared_variant_header(utils::variant_header_deep_copy(header_ptr.get()));

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
    return ITERATOR{m_variant_files, m_variant_headers};
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
   * @brief returns a combined header for the files being read
   * @return a VariantHeader object constructed by combining the headers of the incmoing files
   */
  const inline VariantHeader combined_header() const { return VariantHeader {m_combined_header}; }

 private:
  std::vector<std::shared_ptr<htsFile>> m_variant_files;        ///< vector of the internal file structures of the variant files
  std::vector<std::shared_ptr<bcf_hdr_t>> m_variant_headers;    ///< vector of the internal header structures of the variant files
  std::shared_ptr<bcf_hdr_t> m_combined_header;                 ///< the combined header created from the variant headers

  ///< confirms that the chromosomes in the headers of all of the input files are identical
  // TODO? only handles chromosome names, not lengths
  void validate_header(const std::shared_ptr<bcf_hdr_t>& other_header_ptr) {
    const auto& other_header = VariantHeader{other_header_ptr};
    if (combined_header().chromosomes() != other_header.chromosomes())
      throw HeaderCompatibilityException{"chromosomes in header files are inconsistent"};
  }
};

}  // end namespace gamgee

#endif /* defined(gamgee__multiple_variant_reader__guard) */
