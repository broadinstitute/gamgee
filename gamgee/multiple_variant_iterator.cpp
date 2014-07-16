#include "multiple_variant_iterator.h"

namespace gamgee {

MultipleVariantIterator::MultipleVariantIterator() :
  m_queue {},
  m_variant_vector {}
{}

MultipleVariantIterator::MultipleVariantIterator(const std::vector<vcfFile*> variant_files, const std::shared_ptr<bcf_hdr_t> variant_header) :
  m_queue {},
  m_variant_vector {}
{
  m_variant_vector.reserve(variant_files.size());
  for (auto i = 0u; i < variant_files.size(); i++) {
    m_queue.push(std::make_shared<VariantIterator>(variant_files[i], variant_header));
  }
  fetch_next_vector();
}

MultipleVariantIterator::MultipleVariantIterator(MultipleVariantIterator&& original) :
  m_queue {std::move(original.m_queue)},
  m_variant_vector {std::move(original.m_variant_vector)}
{}

std::vector<Variant>& MultipleVariantIterator::operator*() {
  return m_variant_vector;
}

std::vector<Variant>& MultipleVariantIterator::operator++() {
  fetch_next_vector();
  return m_variant_vector;
}

// NOTE: this method does the minimal work necessary to determine that we have reached the end of iteration
// it is NOT a valid general-purpose inequality method
bool MultipleVariantIterator::operator!=(const MultipleVariantIterator& rhs) {
  return !(m_variant_vector.empty() && rhs.m_variant_vector.empty());
}

bool MultipleVariantIterator::Comparator::operator()(const std::shared_ptr<VariantIterator>& left, const std::shared_ptr<VariantIterator>& right) {
  if ((**left).chromosome() > (**right).chromosome())
    return true;

  if ((**left).chromosome() < (**right).chromosome())
    return false;

  return (**left).alignment_start() > (**right).alignment_start();
}

void MultipleVariantIterator::fetch_next_vector() {
  m_variant_vector.clear();
  auto current_chrom = 0u;
  auto current_start = 0u;

  while (!m_queue.empty()) {
    const auto top_iterator = m_queue.top();
    const auto& variant = **top_iterator;

    if (!m_variant_vector.empty() && !(variant.chromosome() == current_chrom && variant.alignment_start() == current_start))
      break;
    else {
      current_chrom = variant.chromosome();
      current_start = variant.alignment_start();
      m_variant_vector.emplace_back(std::move(variant));
      m_queue.pop();
      top_iterator->operator++();
      if (! top_iterator->empty())
        m_queue.push(std::move(top_iterator));
    }
  }
}

}

