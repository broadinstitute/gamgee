#include "sam_writer.h"

namespace foghorn {

SamWriter::SamWriter(const std::string& output_fname) :
  m_out_file {open_file(output_fname)}
{}

SamWriter::SamWriter(const Sam& sam_record, const std::string& output_fname) :
  m_out_file {open_file(output_fname)},
  m_header{sam_record.header()}
{}

void SamWriter::add_header(const Sam& sam_record) { 
  m_header = SamHeader{sam_record.header()};
  sam_hdr_write(m_out_file, m_header.m_header); 
}

SamWriter::~SamWriter() {
  sam_close(m_out_file);
}

void SamWriter::add_record(const Sam& sam_record) { 
  sam_write1(m_out_file, m_header.m_header, sam_record.body().m_body); 
}

htsFile* SamWriter::open_file(const std::string& output_fname) {
  return hts_open(output_fname.empty() ? "-" : output_fname.c_str(), "w");
}

}
