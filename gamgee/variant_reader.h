#ifndef gamgee__variant_reader__guard
#define gamgee__variant_reader__guard

#include "variant_header.h"
#include "variant_iterator.h"
#include "exceptions.h"

#include "utils/hts_memory.h"
#include "utils/variant_utils.h"

#include "htslib/vcf.h"

#include <string>
#include <fstream>
#include <algorithm>
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
  explicit VariantReader(const std::string& filename) :
    m_variant_file_ptr {},
    m_variant_header_ptr {}
  {
    init_reader(filename);
  }

  /**
   * @brief reads through all records in a file (vcf or bcf) parsing them into Variant
   * objects
   *
   * @param filenames a vector containing a single element: the name of the variant file
   */
  explicit VariantReader(const std::vector<std::string>& filenames) :
    m_variant_file_ptr {},
    m_variant_header_ptr {}
  {
    if (filenames.size() > 1)
      throw SingleInputException{"filenames", filenames.size()};
    if (!filenames.empty())
      init_reader(filenames.front());
  }

  /**
   * @brief reads through all records in a file (vcf or bcf) parsing them into Variant
   * objects but only including the selected samples. To create a sites only file, simply
   * pass an empty vector of samples.
   *
   * @param filename the name of the variant file
   * @param samples the list of samples you want included/excluded from your iteration
   * @param include whether you want these samples to be included or excluded from your iteration.  default = true (include)
   */
  VariantReader(const std::string& filename, const std::vector<std::string>& samples, const bool include = true) :
    m_variant_file_ptr {},
    m_variant_header_ptr {}
  {
    init_reader(filename);
    subset_variant_samples(m_variant_header_ptr.get(), samples, include);
  }

  /**
   * @brief reads through all records in a file (vcf or bcf) parsing them into Variant
   * objects but only including the selected samples. To create a sites only file, simply
   * pass an empty vector of samples.
   *
   * @param filenames a vector containing a single element: the name of the variant file
   * @param samples the list of samples you want included/excluded from your iteration
   * @param include whether you want these samples to be included or excluded from your iteration.  default = true (include)
   */
  VariantReader(const std::vector<std::string>& filenames, const std::vector<std::string>& samples, const bool include = true) :
    m_variant_file_ptr {},
    m_variant_header_ptr {}
  {
    if (filenames.size() > 1)
      throw SingleInputException{"filenames", filenames.size()};
    if (!filenames.empty())
      init_reader(filenames.front());
  }

  /**
   * @brief a VariantReader cannot be copied safely, as it is iterating over a stream.
   */

  VariantReader(const VariantReader& other) = delete;
  VariantReader& operator=(const VariantReader& other) = delete;

  /**
   * @brief a VariantReader can be moved
   */

  VariantReader(VariantReader&& other) = default;
  VariantReader& operator=(VariantReader&& other) = default;

  /**
   * @brief creates a ITERATOR pointing at the start of the input stream (needed by for-each
   * loop)
   *
   * @return a ITERATOR ready to start parsing the file
   */
  ITERATOR begin() const {
    return ITERATOR{ m_variant_file_ptr, m_variant_header_ptr };
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
  std::shared_ptr<htsFile> m_variant_file_ptr;          ///< pointer to the internal file structure of the variant/bam/cram file
  std::shared_ptr<bcf_hdr_t> m_variant_header_ptr;      ///< pointer to the internal header structure of the variant/bam/cram file

  /**
   * @brief initialize the VariantReader (helper function for constructors)
   *
   * @param filename the name of the variant file
   */
  void init_reader (const std::string& filename) {
    auto* file_ptr = bcf_open(filename.empty() ? "-" : filename.c_str(), "r");
    m_variant_file_ptr = utils::make_shared_hts_file (file_ptr);
    m_variant_header_ptr = utils::make_shared_variant_header (bcf_hdr_read (file_ptr));
  }
};

using SingleVariantReader = VariantReader<VariantIterator>;

}  // end of namespace

#endif /* defined(gamgee__variant_reader__guard) */
