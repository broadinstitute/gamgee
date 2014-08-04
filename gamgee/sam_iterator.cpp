#include "sam_iterator.h"
#include "sam.h"
#include "utils/hts_memory.h"

using namespace std;

namespace gamgee {

SamIterator::SamIterator() :
  m_sam_file_ptr {nullptr},
  m_sam_header_ptr {nullptr},
  m_sam_record_ptr {nullptr}
{}

SamIterator::SamIterator(samFile * sam_file_ptr, const std::shared_ptr<bam_hdr_t>& sam_header_ptr) :
  m_sam_file_ptr {sam_file_ptr},
  m_sam_header_ptr {sam_header_ptr},
  m_sam_record_ptr {utils::make_shared_sam(bam_init1())},      ///< important to initialize the record buffer in the constructor so we can reuse it across the iterator
  m_sam_record {fetch_next_record()}
{}

Sam& SamIterator::operator*() {
  return m_sam_record;
}

Sam& SamIterator::operator++() {
  m_sam_record = fetch_next_record();
  return m_sam_record;
}

bool SamIterator::operator!=(const SamIterator& rhs) {
  return m_sam_file_ptr != rhs.m_sam_file_ptr;
}
/**
 * @brief pre-fetches the next sam record
 * @warning we're reusing the existing htslib memory, so users should be aware that all objects from the previous iteration are now stale unless a deep copy has been performed
 */
Sam SamIterator::fetch_next_record() {
 if (sam_read1(m_sam_file_ptr, m_sam_header_ptr.get(), m_sam_record_ptr.get()) < 0) {
    m_sam_file_ptr = nullptr;
    return Sam{};
  }
  return Sam{m_sam_header_ptr, m_sam_record_ptr};
}

}

