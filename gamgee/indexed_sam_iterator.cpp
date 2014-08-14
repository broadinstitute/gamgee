#include "indexed_sam_iterator.h"
#include "sam.h"
#include "utils/hts_memory.h"

#include <iostream> // TODO: delete me post debug

using namespace std;

namespace gamgee {

IndexedSamIterator::IndexedSamIterator() :
  m_sam_file_ptr {nullptr},
  m_sam_index_ptr {nullptr},
  m_sam_header_ptr {nullptr},
  m_interval_list {nullptr},
  m_sam_itr_ptr {nullptr},
  m_sam_record_ptr {nullptr} {
}

IndexedSamIterator::IndexedSamIterator(samFile* sam_file_ptr, hts_idx_t* sam_index_ptr,
    const std::shared_ptr<bam_hdr_t>& sam_header_ptr, const std::vector<std::string>& interval_list) :
  m_sam_file_ptr {sam_file_ptr},
  m_sam_index_ptr {sam_index_ptr},
  m_sam_header_ptr {sam_header_ptr},
  m_interval_list {interval_list},
  m_interval_iterator {m_interval_list.begin()},
  m_sam_itr_ptr {fetch_current_interval()}, // TODO: Handle an empty interval_list. May need to NOT call m_sam_record {fetch_next_record()}?
  m_sam_record_ptr {utils::make_shared_sam(bam_init1())},
  m_sam_record {fetch_next_record()} {
}

IndexedSamIterator::~IndexedSamIterator() {
  if (m_sam_itr_ptr)
    sam_itr_destroy(m_sam_itr_ptr);
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

hts_itr_t* IndexedSamIterator::fetch_current_interval() {
  const auto interval = *m_interval_iterator;
  // For now, use hts_itr_querys to translate gamgee::utils::contig_values.
  return sam_itr_querys(m_sam_index_ptr, m_sam_header_ptr.get(), interval.c_str());
}

Sam IndexedSamIterator::fetch_next_record() {
  std::cout << "m_sam_file_ptr: " << (void*)m_sam_file_ptr << std::endl;
  std::cout << "m_sam_itr_ptr: " << (void*)m_sam_itr_ptr << std::endl;
  std::cout << "last interval: " << *m_interval_iterator << std::endl;
  while (sam_itr_next(m_sam_file_ptr, m_sam_itr_ptr, m_sam_record_ptr.get()) < 0) {
    std::cout << "iterating intervals" << std::endl;
    ++m_interval_iterator;
    if (m_interval_list.end() == m_interval_iterator) {
      std::cout << "exiting intervals" << std::endl;
      m_sam_file_ptr = nullptr;
      return Sam{};
    }
    m_sam_itr_ptr = fetch_current_interval();
    std::cout << "next interval: " << *m_interval_iterator << std::endl;
  }
  std::cout << "returning a Sam" << std::endl;
  std::cout << "m_sam_record_ptr: " << (void*)m_sam_record_ptr.get() << std::endl;
  return Sam{m_sam_header_ptr, m_sam_record_ptr};
}

}
