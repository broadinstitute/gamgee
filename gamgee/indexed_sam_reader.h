#ifndef gamgee__indexed_sam_reader__guard
#define gamgee__indexed_sam_reader__guard

#include "indexed_sam_iterator.h"
#include "utils/hts_memory.h"
#include "utils/sam_utils.h"

#include "htslib/sam.h"

#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>

namespace gamgee {

template<class ITERATOR>

class IndexedSamReader {
  public:

    IndexedSamReader(const std::string& filename, const std::vector<std::string>& interval_list) :
      m_sam_file_ptr {sam_open(filename.c_str(), "r")}, // TODO: Not checking for errors.
      m_sam_index_ptr {sam_index_load(m_sam_file_ptr, filename.c_str())}, // TODO: Not checking for errors.
      m_interval_list {interval_list},
      m_sam_header_ptr { utils::make_shared_sam_header(sam_hdr_read(m_sam_file_ptr)) } // TODO: Not checking for errors.
    {
      std::cout << "In constructor: IndexedSamReader" << std::endl;
      std::cout << "m_sam_file_ptr: " << (void*)m_sam_file_ptr << std::endl;
      std::cout << "m_sam_index_ptr: " << (void*)m_sam_index_ptr << std::endl;
      std::cout << "m_sam_header_ptr: " << (void*)(m_sam_header_ptr.get()) << std::endl;
      std::cout << "Done constructor: IndexedSamReader" << std::endl;
    }

    IndexedSamReader(IndexedSamReader&& other) :
      m_sam_file_ptr {std::move(other.m_sam_file_ptr)},
      m_sam_index_ptr {std::move(other.m_sam_index_ptr)},
      m_interval_list {std::move(other.m_interval_list)},
      m_sam_header_ptr {std::move(other.m_sam_header_ptr)} {
      other.m_sam_file_ptr = nullptr;
    }

    IndexedSamReader& operator=(IndexedSamReader&& other) {
      m_sam_file_ptr = std::move(other.m_sam_file_ptr);
      m_sam_header_ptr = std::move(other.m_sam_header_ptr);
      m_interval_list = std::move(other.m_interval_list);
      other.m_sam_file_ptr = nullptr;
      return *this;
    }

    IndexedSamReader(IndexedSamReader& other) = delete;
    IndexedSamReader& operator=(IndexedSamReader& other) = delete;

    ~IndexedSamReader() {
      if (m_sam_file_ptr) 
        sam_close(m_sam_file_ptr);
      if (m_sam_index_ptr)
        hts_idx_destroy(m_sam_index_ptr); // TODO: Copied from samtools. Is there a {bam,sam}_idx_destroy alias? Didn't find under "BAM/CRAM indexing" in sam.h
    }

    ITERATOR begin() {
      return ITERATOR{m_sam_file_ptr, m_sam_index_ptr, m_sam_header_ptr, m_interval_list};
    }

    ITERATOR end() {
      return ITERATOR{};
    }

    inline SamHeader header() { return SamHeader{m_sam_header_ptr}; }

  private:
    samFile* m_sam_file_ptr;
    hts_idx_t* m_sam_index_ptr;
    std::vector<std::string> m_interval_list;
    std::shared_ptr<bam_hdr_t> m_sam_header_ptr;
};

using IndexedSingleSamReader = IndexedSamReader<IndexedSamIterator>;

}

#endif
