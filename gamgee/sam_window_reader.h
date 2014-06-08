#ifndef __gamgee__sam_window_reader__
#define __gamgee__sam_window_reader__

#include "sam_window_iterator.h"
#include "hts_memory.h"

#include "htslib/sam.h"

#include <string>
#include <iostream>
#include <fstream>
#include <memory>


namespace gamgee {

template<class ITERATOR>
class SamWindowReader {
  public:

    SamWindowReader(const std::string& filename, const uint32_t window_size, const uint32_t step_size) :
      m_sam_file_ptr {sam_open(filename.empty() ? "-" : filename.c_str(), "r")},
      m_sam_header_ptr { make_shared_bam_header(sam_hdr_read(m_sam_file_ptr)) }
    {
      m_window_size = window_size;
      m_step_size = step_size;
    }

    SamWindowReader(SamWindowReader&& other) :
      m_sam_file_ptr {std::move(other.m_sam_file_ptr)},
      m_sam_header_ptr {std::move(other.m_sam_header_ptr)}
    {
      m_window_size = other.m_window_size;
      m_step_size = other.m_step_size;
    }

    ~SamWindowReader() {
      sam_close(m_sam_file_ptr);
    }

    SamWindowReader(const SamWindowReader&) = delete;

    ITERATOR begin() {
      return ITERATOR{m_sam_file_ptr, m_sam_header_ptr, m_window_size, m_step_size};
    }

    ITERATOR end() {
      return ITERATOR{};
    }

    inline SamHeader header() {
      return SamHeader{m_sam_header_ptr};
    }

  private:
    samFile* m_sam_file_ptr;
    const std::shared_ptr<bam_hdr_t> m_sam_header_ptr;
    uint32_t m_window_size;
    uint32_t m_step_size;
};

using SingleSamWindowReader = SamWindowReader<SamWindowIterator>;

}

#endif
