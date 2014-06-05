#include "sam_window_iterator.h"
#include "sam.h"

using namespace std;

namespace gamgee {

SamWindowIterator::SamWindowIterator() :
  m_sam_window{SamWindow{}},
  m_sam_file_ptr{nullptr},
  m_sam_header_ptr{nullptr},
  m_sam_record_ptr{nullptr}
{
  m_is_done_records = true;
  m_current_position = 0;
  m_stop_position = 0;
  m_current_contig = 0;
  m_step_size = 0;
  m_window_size = 0;
}

SamWindowIterator::SamWindowIterator(samFile * sam_file_ptr, const shared_ptr<bam_hdr_t>& sam_header_ptr,
    const uint32_t window_size, const uint32_t step_size) :
  m_sam_window{SamWindow{sam_header_ptr}},
  m_sam_file_ptr{sam_file_ptr},
  m_sam_header_ptr{sam_header_ptr},
  m_sam_record_ptr{make_shared_bam(bam_init1())}
{
  m_is_done_records = false;
  m_current_position = 0;
  m_stop_position = 0;
  m_current_contig = 0;
  m_step_size = step_size;
  m_window_size = window_size;
  move_to_first_window();
}

SamWindowIterator::~SamWindowIterator() {
  m_sam_file_ptr = nullptr;
}

bool SamWindowIterator::is_done_sam_window_iterator() const {
  if (nullptr == m_sam_header_ptr) {
    return true;
  }
  if (m_current_contig >= m_sam_header_ptr.get()->n_targets) {
    return true;
  }
  if (m_current_position >= m_sam_header_ptr.get()->target_len[m_current_contig]) {
    return true;
  }
  return false;
}

SamWindow& SamWindowIterator::operator*() {
  return m_sam_window;
}

SamWindow& SamWindowIterator::operator++() {
  move_to_next_window();
  return m_sam_window;
}

bool SamWindowIterator::operator!=(const SamWindowIterator& rhs) {
  return m_sam_file_ptr != rhs.m_sam_file_ptr &&
      (is_done_sam_window_iterator() != rhs.is_done_sam_window_iterator());
}

void SamWindowIterator::move_to_first_window() {
  if (is_done_sam_window_iterator()) {
    m_is_done_records = true;
    return;
  }
  m_current_position = 0;
  auto target_len = m_sam_header_ptr.get()->target_len[m_current_contig];
  m_stop_position = std::min(target_len, m_window_size);

  retrieve_next_record();
  while (is_keep_appending_records()) {
    append_record();
    retrieve_next_record();
  }
}

void SamWindowIterator::move_to_next_window() {
  auto target_len = m_sam_header_ptr.get()->target_len[m_current_contig];
  auto new_start_position = m_current_position + m_step_size;
  if (new_start_position > target_len) {
    m_current_contig += 1;
    new_start_position = 0;
  }
  m_current_position = new_start_position;
  target_len = m_sam_header_ptr.get()->target_len[m_current_contig];
  m_stop_position = std::min(target_len, m_current_position + m_window_size);
  while (is_keep_dequeuing_records()) {
    dequeue_record();
  }
  while (is_keep_appending_records()) {
    append_record();
    retrieve_next_record();
  }
  m_sam_window.start = m_current_position;
  m_sam_window.stop = m_stop_position;
}

void SamWindowIterator::retrieve_next_record() {
  // WARNING: we're reusing the existing htslib memory, so users should be aware that all
  // objects from the previous iteration are now stale unless a deep copy has been performed
  if (m_is_done_records) {
    return;
  }
  if (sam_read1(m_sam_file_ptr, m_sam_header_ptr.get(), m_sam_record_ptr.get()) < 0) {
    m_is_done_records = true;
    // We're done pulling records. Just finish out the windows.
  }
}

bool SamWindowIterator::is_keep_dequeuing_records() {
  if (m_sam_window.is_empty()) {
    return false;
  }
  auto head = m_sam_window.peek_record().get();
  if (head->core.tid != m_current_contig) {
    return true;
  }
  if (head->core.pos < m_current_position) {
    return true;
  }
  return false;
}

void SamWindowIterator::dequeue_record() {
  m_sam_window.dequeue_record();
}

bool SamWindowIterator::is_keep_appending_records() {
  if (is_done_pulling_records()) {
    return false;
  }
  if (m_sam_record_ptr.get()->core.tid != m_current_contig) {
    return false;
  }
  if (m_sam_record_ptr.get()->core.pos >= m_stop_position) {
    return false;
  }
  return true;
}

void SamWindowIterator::append_record() {
  m_sam_window.enqueue_record(m_sam_record_ptr);
}

bool SamWindowIterator::is_done_pulling_records() {
  return m_is_done_records;
}

}

