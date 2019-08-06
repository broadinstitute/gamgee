#include "sam_writer.h"

#include "../utils/hts_memory.h"

namespace gamgee {

SamWriter::SamWriter(const std::string& output_fname, const bool binary) :
  m_out_file {utils::make_unique_hts_file(open_file(output_fname, binary ? "wb" : "w"))},
  m_header {nullptr}
{}

SamWriter::SamWriter(const SamHeader& header, const std::string& output_fname, const bool binary) :
  m_out_file {utils::make_unique_hts_file(open_file(output_fname, binary ? "wb" : "w"))},
  m_header{header}
{
  write_header();
}

int SamWriter::add_header(const SamHeader& header) {
  m_header = header;
  return write_header();
}

int SamWriter::add_record(const Sam& body) {
  return sam_write1(m_out_file.get(), m_header.m_header.get(), body.m_body.get());
}

htsFile* SamWriter::open_file(const std::string& output_fname, const std::string& mode) {
  return hts_open(output_fname.empty() ? "-" : output_fname.c_str(), mode.c_str());
}

int SamWriter::write_header() const {
  return sam_hdr_write(m_out_file.get(), m_header.m_header.get());
}


}
