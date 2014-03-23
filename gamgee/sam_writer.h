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
  explicit SamWriter(const std::string& output_fname = "-");

  /**
   * @brief Creates a new SamWriter with the header extracted from a Sam record and 
   * using the specified output file name
   *
   * @param sam_record   Sam record to extract the header from
   * @param output_fname file to write to. The default is stdout  (as defined by htslib)
   */
  SamWriter(const Sam& sam_record, const std::string& output_fname = "-");

  void add_record(const Sam& sam_record);
  void add_header(const Sam& sam_record);

  ~SamWriter();

 private:
  htsFile* m_out_file;
  SamHeader m_header;

  static htsFile* open_file(const std::string& output_fname);

};

}
#endif 
