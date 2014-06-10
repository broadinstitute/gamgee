#include "sam_window_node_iterator.h"

namespace gamgee {

SamWindowNodeIterator::SamWindowNodeIterator() :
m_sam_header_ptr{nullptr},
m_sam_record_ptr{nullptr},
m_sam_window_node_ptr{nullptr} {
}

SamWindowNodeIterator::SamWindowNodeIterator(SamWindowNode* sam_window_node_ptr, const std::shared_ptr<bam_hdr_t>& sam_header_ptr) :
m_sam_header_ptr{sam_header_ptr},
m_sam_window_node_ptr{sam_window_node_ptr} {
  m_sam_record = fetch_next_record();
}

bool SamWindowNodeIterator::operator!=(SamWindowNodeIterator& rhs) {
  return m_sam_window_node_ptr != rhs.m_sam_window_node_ptr;
}

Sam& SamWindowNodeIterator::operator*() {
  return m_sam_record;
}

Sam& SamWindowNodeIterator::operator++() {
  m_sam_window_node_ptr = m_sam_window_node_ptr->next;
  m_sam_record = fetch_next_record();
  return m_sam_record;
}

Sam SamWindowNodeIterator::fetch_next_record() {
  if (nullptr == m_sam_window_node_ptr) {
    return Sam{};
  } else {
    m_sam_record_ptr = m_sam_window_node_ptr->sam_record_ptr;
    return Sam{m_sam_record_ptr, m_sam_header_ptr};
  }
}

}
