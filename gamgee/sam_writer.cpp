#include "sam_writer.h"

namespace gamgee {

SamWriter::SamWriter(const std::string& output_fname, const bool binary) :
  m_out_file {open_file(output_fname, binary ? "wb" : "w")},
  m_header {nullptr}
{}

SamWriter::SamWriter(const SamHeader& header, const std::string& output_fname, const bool binary) :
  m_out_file {open_file(output_fname, binary ? "wb" : "w")},
  m_header{header}
{
  write_header();
}

SamWriter::~SamWriter() {
  sam_close(m_out_file);
}

void SamWriter::add_header(const SamHeader& header) { 
  m_header = header;
  write_header();
}

void SamWriter::add_record(const SamBody& body) { 
  sam_write1(m_out_file, m_header.m_header.get(), body.m_body.get());
}

htsFile* SamWriter::open_file(const std::string& output_fname, const std::string& mode) {
  return hts_open(output_fname.empty() ? "-" : output_fname.c_str(), mode.c_str());
}

void SamWriter::write_header() const {
  sam_hdr_write(m_out_file, m_header.m_header.get());
}


}
