#ifndef __gamgee_sam_writer__
#define __gamgee_sam_writer__

#include <string>

#include "sam.h"
#include "sam_header.h"
#include "htslib/sam.h"

namespace gamgee {

class SamWriter {

 public: 
   /**
    * @brief Creates a new SamWriter using the specified output file name
    *
    * @param output_fname file to write to. The default is stdout (as defined by htslib)
    */
  explicit SamWriter(const std::string& output_fname = "-", const bool binary = true);

  /**
   * @brief Creates a new SamWriter with the header extracted from a Sam record and 
   * using the specified output file name
   *
   * @param sam_record   Sam record to extract the header from
   * @param output_fname file to write to. The default is stdout  (as defined by htslib)
   */
  explicit SamWriter(const SamHeader& header, const std::string& output_fname = "-", const bool binary = true);

  void add_record(const SamBody& sam_record);
  void add_header(const SamHeader& sam_record);

  ~SamWriter();

 private:
  htsFile* m_out_file;
  SamHeader m_header;

  static htsFile* open_file(const std::string& output_fname, const std::string& binary);

};

}
#endif 
