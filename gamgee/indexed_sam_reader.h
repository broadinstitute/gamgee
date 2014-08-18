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
      m_sam_header_ptr { utils::make_shared_sam_header(sam_hdr_read(m_sam_file_ptr)) }, // TODO: Not checking for errors.
      m_interval_list {interval_list},
      m_interval_iterator_ptr {nullptr}, // TODO: new? then delete later
      m_sam_itr_ptr_ptr {nullptr} {
        std::cout << "in constructor-a" << std::endl;
      // TODO: Handle an empty interval_list. May need to NOT call m_sam_record {fetch_next_record()}?
      *m_interval_iterator_ptr = m_interval_list.begin();
      std::cout << "in constructor-b" << std::endl;
      *m_sam_itr_ptr_ptr = sam_itr_querys(m_sam_index_ptr, m_sam_header_ptr.get(), (*m_interval_iterator_ptr)->c_str());
      std::cout << "in constructor" << std::endl;
    }

    IndexedSamReader(IndexedSamReader&& other) :
      m_sam_file_ptr {std::move(other.m_sam_file_ptr)},
      m_sam_index_ptr {std::move(other.m_sam_index_ptr)},
      m_sam_header_ptr {std::move(other.m_sam_header_ptr)},
      m_interval_list {std::move(other.m_interval_list)},
      m_interval_iterator_ptr {std::move(other.m_interval_iterator_ptr)},
      m_sam_itr_ptr_ptr {std::move(other.m_sam_itr_ptr_ptr)} {
      other.m_sam_file_ptr = nullptr;
      other.m_sam_index_ptr = nullptr;
      other.m_sam_itr_ptr_ptr = nullptr;
    }

    IndexedSamReader& operator=(IndexedSamReader&& other) {
      m_sam_file_ptr = std::move(other.m_sam_file_ptr);
      m_sam_index_ptr = std::move(other.m_sam_index_ptr);
      m_sam_header_ptr = std::move(other.m_sam_header_ptr);
      m_interval_list = std::move(other.m_interval_list);
      m_interval_iterator_ptr = std::move(other.m_interval_iterator_ptr);
      m_sam_itr_ptr_ptr = std::move(other.m_sam_itr_ptr_ptr);
      other.m_sam_file_ptr = nullptr;
      other.m_sam_index_ptr = nullptr;
      other.m_sam_itr_ptr_ptr = nullptr;
      return *this;
    }

    IndexedSamReader(IndexedSamReader& other) = delete;
    IndexedSamReader& operator=(IndexedSamReader& other) = delete;

    ~IndexedSamReader() {
      if (*m_sam_itr_ptr_ptr)
        sam_itr_destroy(*m_sam_itr_ptr_ptr);
      if (m_sam_index_ptr)
        hts_idx_destroy(m_sam_index_ptr); // TODO: Copied from samtools. Is there a {bam,sam}_idx_destroy alias? Didn't find under "BAM/CRAM indexing" in sam.h
      if (m_sam_file_ptr) 
        sam_close(m_sam_file_ptr);
    }

    ITERATOR begin() {
      return ITERATOR{m_sam_file_ptr, m_sam_index_ptr, m_sam_header_ptr, m_interval_list, m_interval_iterator_ptr, m_sam_itr_ptr_ptr};
    }

    ITERATOR end() {
      return ITERATOR{};
    }

    inline SamHeader header() { return SamHeader{m_sam_header_ptr}; }

  private:
    samFile* m_sam_file_ptr;
    hts_idx_t* m_sam_index_ptr;
    std::shared_ptr<bam_hdr_t> m_sam_header_ptr;
    std::vector<std::string> m_interval_list;
    std::vector<std::string>::iterator* m_interval_iterator_ptr;
    hts_itr_t ** m_sam_itr_ptr_ptr; // TODO: Shared with the iterator... so should this be a shared_ptr or is that not necessary?
};

using IndexedSingleSamReader = IndexedSamReader<IndexedSamIterator>;

}

#endif
