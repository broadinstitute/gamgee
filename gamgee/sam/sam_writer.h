#ifndef gamgee__sam_writer__guard
#define gamgee__sam_writer__guard

#include <string>
#include <memory>

#include "sam.h"
#include "sam_header.h"

#include "../utils/hts_memory.h"

#include "htslib/sam.h"

namespace gamgee {

/**
 * @brief utility class to write out a SAM/BAM/CRAM file to any stream
 * @todo add serialization option
 */
class SamWriter {

 public: 

  /**
   * @brief Creates a new SamWriter using the specified output file name
   * @param output_fname file to write to. The default is stdout (as defined by htslib)
   * @param binary whether the output should be in BAM (true) or SAM format (false) 
   * @note the header is copied and managed internally
   */
  explicit SamWriter(const std::string& output_fname = "-", const bool binary = true);

  /**
   * @brief Creates a new SamWriter with the header extracted from a Sam record and using the specified output file name
   * @param header       SamHeader object to make a copy from
   * @param output_fname file to write to. The default is stdout  (as defined by htslib)
   * @param binary whether the output should be in BAM (true) or SAM format (false) 
   * @note the header is copied and managed internally
   */
  explicit SamWriter(const SamHeader& header, const std::string& output_fname = "-", const bool binary = true);

  /**
   * @brief a SamWriter cannot be copied safely, as it is iterating over a stream.
   */

  SamWriter(const SamWriter& other) = delete;
  SamWriter& operator=(const SamWriter& other) = delete;

  /**
   * @brief a SamWriter can be moved
   */

  SamWriter(SamWriter&& other) = default;
  SamWriter& operator=(SamWriter&& other) = default;

  /**
   * @brief Adds a record to the file stream
   * @param body the record
   */
  int add_record(const Sam& body);

  /**
   * @brief Adds a header to the file stream.
   * @param header the header
   * @note the header is a requirement to add records
   */
  int add_header(const SamHeader& header);

 private:
  std::unique_ptr<htsFile, utils::HtsFileDeleter> m_out_file;  ///< the file or stream to write out to ("-" means stdout)
  SamHeader m_header;                   ///< holds a copy of the header throughout the production of the output (necessary for every record that gets added)

  static htsFile* open_file(const std::string& output_fname, const std::string& binary);
  int write_header() const;

};

}

#endif // gamgee__sam_writer__guard
