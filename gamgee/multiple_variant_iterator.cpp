#include "multiple_variant_iterator.h"

namespace gamgee {

MultipleVariantIterator::MultipleVariantIterator() :
  m_queue {},
  m_variant_vector {}
{}

// TODO: passthrough
MultipleVariantIterator::MultipleVariantIterator(const std::vector<vcfFile*> variant_files, const std::vector<const std::shared_ptr<bcf_hdr_t> > variant_headers) :
  m_queue {},
  m_variant_vector {}
{
  m_queue.push(std::unique_ptr<VariantIterator>{new VariantIterator{variant_files[0], variant_headers[0]}});
  m_variant_vector = fetch_next_record();
}

MultipleVariantIterator::MultipleVariantIterator(MultipleVariantIterator&& original) :
  m_queue {std::move(original.m_queue)},
  m_variant_vector {std::move(original.m_variant_vector)}
{}

std::vector<Variant>& MultipleVariantIterator::operator*() {
  return m_variant_vector;
}

std::vector<Variant>& MultipleVariantIterator::operator++() {
  m_variant_vector = fetch_next_record();
  return m_variant_vector;
}

bool MultipleVariantIterator::operator!=(const MultipleVariantIterator& rhs) {
  // TODO: I can't compare the queues due to const-ness, so I have to use a pointer instead?  Is that normal?
  return &m_queue != &rhs.m_queue;
}

/**
 * @brief pre-fetches the next variant record
 */
// TODO: temp passthrough
std::vector<Variant> MultipleVariantIterator::fetch_next_record() {
  if (m_queue.empty())
    return std::vector<Variant>{};
  else {
    auto nextIter = m_queue.top();
    auto retVal = std::vector<Variant>{**nextIter};
    m_queue.pop();
    auto newIter = nextIter;	// de-const
    newIter->operator++();
    m_queue.push(std::move(newIter));
    return retVal;
  }
}

}

