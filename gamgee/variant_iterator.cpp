#include "variant_iterator.h"
#include "variant.h"
#include "utils/hts_memory.h"

using namespace std;

namespace gamgee {

VariantIterator::VariantIterator() :
  m_variant_file_ptr {nullptr},
  m_variant_header_ptr {nullptr},
  m_variant_record_ptr {nullptr}
{}

VariantIterator::VariantIterator(vcfFile* variant_file_ptr, const std::shared_ptr<bcf_hdr_t>& variant_header_ptr) :
  m_variant_file_ptr {variant_file_ptr},
  m_variant_header_ptr {variant_header_ptr},
  m_variant_record_ptr {utils::make_shared_variant(bcf_init1())},      ///< important to initialize the record buffer in the constructor so we can reuse it across the iterator
  m_variant_record {fetch_next_record()}
{}

VariantIterator::VariantIterator(VariantIterator&& original) :
  m_variant_file_ptr   {move(original.m_variant_file_ptr)},
  m_variant_header_ptr {move(original.m_variant_header_ptr)},
  m_variant_record_ptr {move(original.m_variant_record_ptr)},
  m_variant_record     {move(original.m_variant_record)}
{}

Variant& VariantIterator::operator*() {
  return m_variant_record;
}

Variant& VariantIterator::operator++() {
  m_variant_record = fetch_next_record();
  return m_variant_record;
}

bool VariantIterator::operator!=(const VariantIterator& rhs) {
  return m_variant_file_ptr != rhs.m_variant_file_ptr;
}
/**
 * @brief pre-fetches the next variant record
 * @warning we're reusing the existing htslib memory, so users should be aware that all objects from the previous iteration are now stale unless a deep copy has been performed
 */
Variant VariantIterator::fetch_next_record() {
 if (bcf_read1(m_variant_file_ptr, m_variant_header_ptr.get(), m_variant_record_ptr.get()) < 0) {
    m_variant_file_ptr = nullptr;
    return Variant{};
  }
  return Variant{m_variant_header_ptr, m_variant_record_ptr};
}

}

