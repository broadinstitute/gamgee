#ifndef gamgee__reference_iterator__guard
#define gamgee__reference_iterator__guard

#include "fastq_iterator.h"
#include "fastq_reader.h"

namespace gamgee {

/**
 * @brief Utility class to access reference bases in a FastA-formatted reference genome
 *
 * @warn the chromosomes in this reference must be accessed in order, or an exception will be thrown
 */
class ReferenceIterator {
 public:
  ReferenceIterator(const std::string& filename) :
    m_iterator {FastqReader{filename}.begin()},
    m_sequence {*m_iterator}
  {}

  /**
   * @brief return the reference base character at the desired location
   * @param chromosome the chromosome of the desired base
   * @param one_based_location the one-based genomic location of the base
   */
  const char ref_base(const std::string& chromosome, const int one_based_location);

 private:
  FastqIterator m_iterator;     ///< @brief the current state of the iterator through the FastA input file
  Fastq& m_sequence;            ///< @brief a reference to the iterator's current sequence
};

} // namespace gamgee

#endif /* gamgee__reference_iterator__guard */
