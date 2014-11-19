#include "multiple_variant_iterator.h"

namespace gamgee {

MultipleVariantIterator::MultipleVariantIterator(const std::vector<std::shared_ptr<htsFile>>& variant_files, const std::vector<std::shared_ptr<bcf_hdr_t>>& variant_headers) :
  m_queue {},
  m_variant_vector {}
{
  m_variant_vector.reserve(variant_files.size());
  for (auto i = 0u; i < variant_files.size(); i++) {
    m_queue.push(VariantIteratorIndexPair{std::make_shared<VariantIterator>(variant_files[i], variant_headers[i]), i});
  }
  fetch_next_vector();
}

std::vector<VariantIndexPair>& MultipleVariantIterator::operator*() {
  return m_variant_vector;
}

std::vector<VariantIndexPair>& MultipleVariantIterator::operator++() {
  fetch_next_vector();
  return m_variant_vector;
}

// NOTE: this method does the minimal work necessary to determine that we have reached the end of iteration
// it is NOT a valid general-purpose inequality method
bool MultipleVariantIterator::operator!=(const MultipleVariantIterator& rhs) {
  return !(m_variant_vector.empty() && rhs.m_variant_vector.empty());
}

bool MultipleVariantIterator::Comparator::operator()(const VariantIteratorIndexPair& left, const VariantIteratorIndexPair& right) {
  if ((**(left.first)).chromosome() > (**(right.first)).chromosome())
    return true;

  if ((**(left.first)).chromosome() < (**(right.first)).chromosome())
    return false;

  return (**(left.first)).alignment_start() > (**(right.first)).alignment_start();
}

void MultipleVariantIterator::fetch_next_vector() {
  m_variant_vector.clear();
  auto current_chrom = 0u;
  auto current_start = 0u;

  while (!m_queue.empty()) {
    const auto top_queue_elem = m_queue.top();
    const auto top_iterator = top_queue_elem.first;
    const auto& variant = **top_iterator;

    if (!m_variant_vector.empty() && !(variant.chromosome() == current_chrom && variant.alignment_start() == current_start))
      break;
    else {
      current_chrom = variant.chromosome();
      current_start = variant.alignment_start();
      m_variant_vector.push_back(VariantIndexPair{variant, top_queue_elem.second});

      m_queue.pop();
      top_iterator->operator++();
      if (! top_iterator->empty())
        m_queue.push(std::move(top_queue_elem));
    }
  }
}

}

