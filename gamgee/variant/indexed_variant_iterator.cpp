#include "indexed_variant_iterator.h"
#include "variant_iterator.h"

#include "htslib/vcf.h"

#include <memory>
#include <string>
#include <vector>

using namespace std;

namespace gamgee {
const std::vector<std::string> IndexedVariantIterator::all_intervals = {"."};

IndexedVariantIterator::IndexedVariantIterator() :
  VariantIterator {},
  m_variant_index_ptr {},
  m_interval_list {},
  m_interval_iter {},
  m_index_iter_ptr {}
  {}

IndexedVariantIterator::IndexedVariantIterator(const std::shared_ptr<htsFile>& file_ptr,
                                               const std::shared_ptr<hts_idx_t>& index_ptr,
                                               const std::shared_ptr<bcf_hdr_t>& header_ptr,
                                               const std::vector<std::string>& interval_list) :
  VariantIterator { file_ptr, header_ptr },
  m_variant_index_ptr { index_ptr },
  m_interval_list { interval_list.empty() ? all_intervals : interval_list },
  m_interval_iter { m_interval_list.begin() },
  m_index_iter_ptr { utils::make_unique_hts_itr(bcf_itr_querys(m_variant_index_ptr.get(), m_variant_header_ptr.get(), m_interval_iter->c_str())) }
{
  fetch_next_record();
}

bool IndexedVariantIterator::operator!=(const IndexedVariantIterator& rhs) {
  return m_variant_file_ptr != rhs.m_variant_file_ptr &&
    m_index_iter_ptr != rhs.m_index_iter_ptr;
}

/**
 * @brief pre-fetches the next variant record
 * @warning we're reusing the existing htslib memory, so users should be aware that all objects from the previous iteration are now stale unless a deep copy has been performed
 */
void IndexedVariantIterator::fetch_next_record() {
  while (bcf_itr_next(m_variant_file_ptr, m_index_iter_ptr.get(), m_variant_record_ptr.get()) < 0) {
    ++m_interval_iter;
    if (m_interval_list.end() == m_interval_iter) {
      m_variant_file_ptr.reset();
      m_variant_record = Variant{};
      return;
    }
    m_index_iter_ptr.reset(bcf_itr_querys(m_variant_index_ptr.get(), m_variant_header_ptr.get(), m_interval_iter->c_str()));
  }
}

}

