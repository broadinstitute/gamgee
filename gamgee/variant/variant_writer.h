#ifndef gamgee__variant_writer__guard
#define gamgee__variant_writer__guard

#include <string>
#include <memory>
#include <zlib.h>

#include "variant.h"
#include "variant_header.h"

#include "../utils/hts_memory.h"

#include "htslib/vcf.h"

namespace gamgee {

/**
 * @brief utility class to write out a VCF/BCF file to any stream
 * @todo add serialization option
 */
class VariantWriter {

 public: 

  /**
   * @brief Creates a new VariantWriter using the specified output file name
   * @param output_fname file to write to. The default is stdout (as defined by htslib)
   * @param binary whether the output should be in BCF (true) or VCF format (false)
   * @param compression_level optional zlib compression level. 0 for none, 1 for best speed, 9 for best compression
   * @note the header is copied and managed internally
   */
  explicit VariantWriter(const std::string& output_fname = "-", const bool binary = true, const int compression_level = Z_DEFAULT_COMPRESSION);

  /**
   * @brief Creates a new VariantWriter with the header extracted from a Variant record and using the specified output file name
   * @param header a VariantHeader object to make a copy from
   * @param output_fname file to write to. The default is stdout  (as defined by htslib)
   * @param binary whether the output should be in BCF (true) or VCF format (false)
   * @param compression_level optional zlib compression level. 0 for none, 1 for best speed, 9 for best compression
   * @note the header is copied and managed internally
   */
  explicit VariantWriter(const VariantHeader& header, const std::string& output_fname = "-", const bool binary = true, const int compression_level = Z_DEFAULT_COMPRESSION);

  /**
   * @brief a VariantWriter cannot be copied safely, as it is iterating over a stream.
   */

  VariantWriter(const VariantWriter& other) = delete;
  VariantWriter& operator=(const VariantWriter& other) = delete;

  /**
   * @brief a VariantWriter can be moved
   */

  VariantWriter(VariantWriter&& other) = default;
  VariantWriter& operator=(VariantWriter&& other) = default;

  /**
   * @brief Adds a record to the file stream
   * @param body the record
   */
  void add_record(const Variant& body);

  /**
   * @brief Adds a header to the file stream.
   * @param header the header
   * @note the header is a requirement to add records
   */
  void add_header(const VariantHeader& header);

 private:
  std::unique_ptr<htsFile, utils::HtsFileDeleter> m_out_file;  ///< the file or stream to write out to ("-" means stdout)
  VariantHeader m_header;               ///< holds a copy of the header throughout the production of the output (necessary for every record that gets added)

  static htsFile* open_file(const std::string& output_fname, const std::string& binary);
  void write_header() const;
  std::string write_mode(const bool binary, const int compression_level) const;
};

}

#endif // gamgee__variant_writer__guard
