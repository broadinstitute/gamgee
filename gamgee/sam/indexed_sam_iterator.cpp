#include "indexed_sam_iterator.h"
#include "sam.h"

#include "../utils/hts_memory.h"

using namespace std;

namespace gamgee {

IndexedSamIterator::IndexedSamIterator() :
  m_sam_file_ptr {nullptr},
  m_sam_index_ptr {nullptr},
  m_sam_header_ptr {nullptr},
  m_sam_itr_ptr {nullptr},
  m_sam_record_ptr {nullptr} {
}

IndexedSamIterator::IndexedSamIterator(const std::shared_ptr<htsFile>& sam_file_ptr, const std::shared_ptr<hts_idx_t>& sam_index_ptr,
    const std::shared_ptr<bam_hdr_t>& sam_header_ptr, const std::vector<std::string>& interval_list) :
  m_sam_file_ptr {sam_file_ptr},
  m_sam_index_ptr {sam_index_ptr},
  m_sam_header_ptr {sam_header_ptr},
  m_interval_list {interval_list},
  m_interval_iterator {m_interval_list.begin()},
  m_sam_itr_ptr {utils::make_unique_hts_itr(sam_itr_querys(m_sam_index_ptr.get(), m_sam_header_ptr.get(), (*m_interval_iterator).c_str()))},
  m_sam_record_ptr {utils::make_shared_sam(bam_init1())},
  m_sam_record {m_sam_header_ptr, m_sam_record_ptr} {
    fetch_next_record();
}

Sam& IndexedSamIterator::operator*() {
  return m_sam_record;
}

Sam& IndexedSamIterator::operator++() {
  fetch_next_record();
  return m_sam_record;
}

bool IndexedSamIterator::operator!=(const IndexedSamIterator& rhs) {
  return m_sam_file_ptr != rhs.m_sam_file_ptr;
}

void IndexedSamIterator::fetch_next_record() {
  while (sam_itr_next(m_sam_file_ptr.get(), m_sam_itr_ptr.get(), m_sam_record_ptr.get()) < 0) {
    ++m_interval_iterator;
    if (m_interval_list.end() == m_interval_iterator) {
      m_sam_file_ptr = nullptr;
      return;
    }
    m_sam_itr_ptr.reset(sam_itr_querys(m_sam_index_ptr.get(), m_sam_header_ptr.get(), (*m_interval_iterator).c_str()));
  }
}

const std::string& IndexedSamIterator::current_interval() const{
  return *m_interval_iterator;
}


}
