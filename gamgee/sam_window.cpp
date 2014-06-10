#include "sam_window.h"

namespace gamgee {

SamWindow::SamWindow() :
  m_sam_header_ptr{nullptr},
  m_head{nullptr},
  m_tail{nullptr},
  m_size{0}
{
}

SamWindow::SamWindow(const std::shared_ptr<bam_hdr_t>& sam_header_ptr) :
  m_sam_header_ptr{sam_header_ptr},
  m_head{nullptr},
  m_tail{nullptr},
  m_size{0}
{
}

SamWindow::~SamWindow() {
  while (!is_empty()) {
    dequeue_record();
  }
}

std::shared_ptr<bam1_t> SamWindow::peek_record() {
  return m_head->sam_record_ptr;
}

std::shared_ptr<bam1_t> SamWindow::dequeue_record() {
  auto to_delete = m_head;
  auto to_return = m_head->sam_record_ptr;
  if (m_head == m_tail) {
    m_tail = nullptr;
  }
  m_head = m_head->next;
  --m_size;
  delete(to_delete);
  return to_return;
}

void SamWindow::enqueue_record(std::shared_ptr<bam1_t> sam_record_ptr) {
  if (m_tail == nullptr) {
    m_head = new SamWindowNode{};
    m_head->sam_record_ptr = std::shared_ptr<bam1_t>(bam_deep_copy(sam_record_ptr.get()));
    m_head->next = nullptr;
    m_tail = m_head;
  } else {
    m_tail->next = new SamWindowNode{};
    m_tail->next->sam_record_ptr = std::shared_ptr<bam1_t>(bam_deep_copy(sam_record_ptr.get()));
    m_tail->next->next = nullptr;
    m_tail = m_tail->next;
  }
  ++m_size;
}

bool SamWindow::is_empty() {
  return m_size == 0;
}

int32_t SamWindow::size() const {
  return m_size;
}

SamWindowNodeIterator SamWindow::begin() const {
  return SamWindowNodeIterator{m_head, m_sam_header_ptr};
}

SamWindowNodeIterator SamWindow::end() const {
  return SamWindowNodeIterator{};
}

}
