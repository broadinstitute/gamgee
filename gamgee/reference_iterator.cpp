#include "reference_iterator.h"

#include "exceptions.h"
#include "fastq_iterator.h"

namespace gamgee {

const char ReferenceIterator::ref_base(const std::string& chromosome, const int one_based_location) {
  while (m_iterator != FastqIterator{} && chromosome != m_sequence.name())
    m_sequence = ++m_iterator;

  if (m_iterator == FastqIterator{})
    throw new ChromosomeNotFoundException{chromosome};

  if (one_based_location > m_sequence.sequence().size())
    throw new ChromosomeSizeException{chromosome, m_sequence.sequence().size(), one_based_location};

  return m_sequence.sequence()[one_based_location-1];
}

} // namespace gamgee

