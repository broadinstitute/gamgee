#include "indexed_sam_iterator.h"
#include "sam.h"
#include "utils/hts_memory.h"

using namespace std;

namespace gamgee {

IndexedSamIterator::IndexedSamIterator() :
  m_sam_file_ptr {nullptr},
  m_sam_index_ptr {nullptr},
  m_sam_header_ptr {nullptr},
  m_sam_itr_ptr_ptr {nullptr},
  m_sam_record_ptr {nullptr} {
}

IndexedSamIterator::IndexedSamIterator(samFile* sam_file_ptr, hts_idx_t* sam_index_ptr, const std::shared_ptr<bam_hdr_t>& sam_header_ptr,
    const std::vector<std::string>& interval_list, std::vector<std::string>::iterator* interval_iterator_ptr,
    hts_itr_t** sam_itr_ptr_ptr) :
  m_sam_file_ptr {sam_file_ptr},
  m_sam_index_ptr {sam_index_ptr},
  m_sam_header_ptr {sam_header_ptr},
  m_interval_list {interval_list},
  m_interval_iterator_ptr {interval_iterator_ptr},
  m_sam_itr_ptr_ptr {sam_itr_ptr_ptr},
  m_sam_record_ptr {utils::make_shared_sam(bam_init1())}, // TODO: Does this need to be released?
  m_sam_record {fetch_next_record()} {
}

IndexedSamIterator::IndexedSamIterator(IndexedSamIterator&& other) :
  m_sam_file_ptr {std::move(other.m_sam_file_ptr)},
  m_sam_index_ptr {std::move(other.m_sam_index_ptr)},
  m_sam_header_ptr {std::move(other.m_sam_header_ptr)},
  m_interval_list {std::move(other.m_interval_list)},
  m_interval_iterator_ptr {std::move(other.m_interval_iterator_ptr)},
  m_sam_itr_ptr_ptr {std::move(other.m_sam_itr_ptr_ptr)},
  m_sam_record_ptr {std::move(other.m_sam_record_ptr)},
  m_sam_record {std::move(other.m_sam_record)} {
  other.m_sam_file_ptr = nullptr;
}

IndexedSamIterator& IndexedSamIterator::operator=(IndexedSamIterator&& other) {
  m_sam_file_ptr = std::move(other.m_sam_file_ptr);
  m_sam_index_ptr = std::move(other.m_sam_index_ptr);
  m_sam_header_ptr = std::move(other.m_sam_header_ptr);
  m_interval_list = std::move(other.m_interval_list);
  m_interval_iterator_ptr = std::move(other.m_interval_iterator_ptr);
  m_sam_itr_ptr_ptr = std::move(other.m_sam_itr_ptr_ptr);
  m_sam_record_ptr = std::move(other.m_sam_record_ptr);
  m_sam_record = std::move(other.m_sam_record);
  other.m_sam_file_ptr = nullptr;
  return *this;
}

Sam& IndexedSamIterator::operator*() {
  return m_sam_record;
}

Sam& IndexedSamIterator::operator++() {
  m_sam_record = fetch_next_record();
  return m_sam_record;
}

bool IndexedSamIterator::operator!=(const IndexedSamIterator& rhs) {
  return m_sam_file_ptr != rhs.m_sam_file_ptr;
}

Sam IndexedSamIterator::fetch_next_record() {
  while (sam_itr_next(m_sam_file_ptr, *m_sam_itr_ptr_ptr, m_sam_record_ptr.get()) < 0) {
    ++(*m_interval_iterator_ptr);
    if (m_interval_list.end() == (*m_interval_iterator_ptr)) {
      m_sam_file_ptr = nullptr;
      return Sam{};
    }
    sam_itr_destroy(*m_sam_itr_ptr_ptr);
    *m_sam_itr_ptr_ptr = sam_itr_querys(m_sam_index_ptr, m_sam_header_ptr.get(), (*m_interval_iterator_ptr)->c_str());
  }
  return Sam{m_sam_header_ptr, m_sam_record_ptr};
}

}
