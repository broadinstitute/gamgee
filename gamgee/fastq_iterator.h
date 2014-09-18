#ifndef gamgee__fastq_iterator__guard
#define gamgee__fastq_iterator__guard

#include "fastq.h"

#include <memory>
#include <istream>

namespace gamgee {

/**
 * @brief Utility class to enable for-each style iteration in the FastqReader class
 */
class FastqIterator {
 public:

  /**
    * @brief creates an empty iterator (used for the end() method) 
    */
  FastqIterator();

  /**
    * @brief initializes a new iterator based on an input stream (e.g. fastq/a file, stdin, ...)
    *
    * @param in input stream
    */
  explicit FastqIterator(std::shared_ptr<std::istream>& in);

  /**
    * @brief a FastqIterator should never be copied as the underlying stream can only be
    * manipulated by one object.
    */
  FastqIterator(const FastqIterator&) = delete;
  FastqIterator& operator=(const FastqIterator&) = delete;

  /**
    * @brief a FastqIterator move constructor guarantees all objects will have the same state.
    */
  FastqIterator(FastqIterator&&) = default;
  FastqIterator& operator=(FastqIterator&&) = default;

  /**
    * @brief equality operator
    *
    * @param rhs the other FastqIterator to compare to
    *
    * @return whether or not the two iterators are the same (e.g. have the same input stream on the same
    * status)
    */
  bool operator==(const FastqIterator& rhs) const;

  /**
    * @brief inequality operator (needed by for-each loop)
    *
    * @param rhs the other FastqIterator to compare to
    *
    * @return whether or not the two iterators are the same (e.g. have the same input stream on the same
    * status)
    */
  bool operator!=(const FastqIterator& rhs) const;

  /**
    * @brief dereference operator (needed by for-each loop)
    *
    * @return a reference to a parsed Fastq object (current in the stream)
    */
  Fastq& operator*();

  /**
    * @brief increment operator (needed by for-each loop)
    *
    * @return the next parsed Fastq object (next in the stream)
    */
  Fastq& operator++();
  
 private:
  std::shared_ptr<std::istream> m_input_stream;         ///< a pointer to the input stream
  Fastq m_element;              ///< the current parsed fastq/fasta element
  bool m_is_fastq;              ///< whether we are parsing fastq's or fasta's from the input stream
  char m_eos_delim;             ///< delimits the end of the sequence field in the fastq/fasta file
  char m_bor_delim;             ///< delimits the beginning of the record in the fastq/fasta file
  
  Fastq fetch_next_element();
  std::string parse_comment();
  std::string parse_seq();
  std::string parse_quals(uint32_t seq_length);
  const std::string parse_multiline();
  void skip_new_lines();
};

}  // end namespace gamgee

#endif // gamgee__fastq_iterator__guard
