#include "variant_writer.h"

#include "utils/hts_memory.h"

namespace gamgee {

VariantWriter::VariantWriter(const std::string& output_fname, const bool binary) :
  m_out_file {utils::make_shared_hts_file(open_file(output_fname, binary ? "wb" : "w"))},
  m_header {nullptr}
{}

VariantWriter::VariantWriter(const VariantHeader& header, const std::string& output_fname, const bool binary) :
  m_out_file {utils::make_shared_hts_file(open_file(output_fname, binary ? "wb" : "w"))},
  m_header{header}
{
  write_header();
}

void VariantWriter::add_header(const VariantHeader& header) { 
  m_header = header;
  write_header();
}

void VariantWriter::add_record(const Variant& body) { 
  bcf_write1(m_out_file.get(), m_header.m_header.get(), body.m_body.get());
}

htsFile* VariantWriter::open_file(const std::string& output_fname, const std::string& mode) {
  return hts_open(output_fname.empty() ? "-" : output_fname.c_str(), mode.c_str());
}

void VariantWriter::write_header() const {
  bcf_hdr_write(m_out_file.get(), m_header.m_header.get());
}


}
