#include "sam_writer.h"

namespace gamgee {

SamWriter::SamWriter(const std::string& output_fname) :
  m_out_file {open_file(output_fname)}
{}

SamWriter::SamWriter(const SamHeader& sam, const std::string& output_fname) :
  m_out_file {open_file(output_fname)},
  m_header{sam}
{}

void SamWriter::add_header(const SamHeader& header) { 
  m_header = SamHeader{header};
  sam_hdr_write(m_out_file, m_header.m_header); 
}

SamWriter::~SamWriter() {
  sam_close(m_out_file);
}

void SamWriter::add_record(const SamBody& body) { 
  sam_write1(m_out_file, m_header.m_header, body.m_body); 
}

htsFile* SamWriter::open_file(const std::string& output_fname) {
  return hts_open(output_fname.empty() ? "-" : output_fname.c_str(), "w");
}

}
