#include "htslib/vcf.h"
#include "htslib/synced_bcf_reader.h"

#include "synced_variant_iterator.h"
#include "variant.h"

using namespace std;

namespace gamgee {

SyncedVariantIterator::SyncedVariantIterator() :
  m_synced_readers {},
  m_variant_vector {}
{}

SyncedVariantIterator::SyncedVariantIterator(const std::shared_ptr<bcf_srs_t>& synced_readers) :
  m_synced_readers {synced_readers},
  m_variant_vector {}
{
  m_variant_vector.reserve(m_synced_readers->nreaders);
  fetch_next_record();
}

std::vector<Variant>& SyncedVariantIterator::operator*() {
  return m_variant_vector;
}

std::vector<Variant>& SyncedVariantIterator::operator++() {
  fetch_next_record();
  return m_variant_vector;
}

// NOTE: this method does the minimal work necessary to determine that we have reached the end of iteration
// it is NOT a valid general-purpose inequality method
bool SyncedVariantIterator::operator!=(const SyncedVariantIterator& rhs) {
  return !(m_variant_vector.empty() && rhs.m_variant_vector.empty());
}

/**
 * @brief initializes m_headers_vector
 */
void SyncedVariantIterator::init_headers_vector() {
  m_headers_vector.reserve(m_synced_readers->nreaders);
  for (int idx = 0; idx < m_synced_readers->nreaders; idx++) {
    auto* hdr_ptr = utils::variant_header_deep_copy(bcf_sr_get_header(m_synced_readers.get(), idx));
    m_headers_vector.emplace_back(utils::make_shared_variant_header(hdr_ptr));
  }
}

/**
 * @brief pre-fetches the next variant record
 */
void SyncedVariantIterator::fetch_next_record() {
  if (!bcf_sr_next_line(m_synced_readers.get())) {
    m_variant_vector.clear();
    return;
  }
  else {
    m_variant_vector.clear();
    m_variant_vector.resize(m_synced_readers->nreaders);
  }
    
  // can't initialize until a line has been read
  if (m_headers_vector.empty())
    init_headers_vector();

  for (int idx = 0; idx < m_synced_readers->nreaders; idx++) {
    if (bcf_sr_has_line(m_synced_readers.get(), idx)) {
      // can't cache variant bodies because they may change location in the synced reader
      auto* body_ptr = utils::variant_deep_copy(bcf_sr_get_line(m_synced_readers.get(), idx));
      m_variant_vector[idx] = Variant{m_headers_vector[idx], utils::make_shared_variant(body_ptr)};
    }
  }
}

}

