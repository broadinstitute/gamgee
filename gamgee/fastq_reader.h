#ifndef gamgee__fastq_reader__guard
#define gamgee__fastq_reader__guard

#include "fastq_iterator.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

namespace gamgee {

/**
 * @brief Utility class to read many Fastq records from a stream (e.g. Fastq file, stdin, ...) in a
 * for-each loop in a for-each loop.
 *
 * This class is designed to parse fastq files in for-each loops with the following signature:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (auto& record : FastqReader(filename))
 *   do_something_with_fastq(record);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * You can also use it with the stdin or any other stream by simply giving the reference to the
 * stream in the constructor, like so:
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (auto& record : FastqReader(&std::cin))
 *   do_something_with_fastq(record);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Although one could use it as an iterator, if your goal is to do so, you should use the FastqIterator
 * class
 */
class FastqReader {
 public:

  /**
    * @brief reads through all records in a file (fasta or fastq) parsing them into Fastq
    * objects
    *
    * @param filename the name of the fasta/fastq file
    */
  explicit FastqReader(const std::string& filename);

  /**
    * @brief reads through all records in a file (fasta or fastq) parsing them into Fastq
    * objects
    *
    * @param filenames a vector containing a single element: the name of the fasta/fastq file
    */
  explicit FastqReader(const std::vector<std::string>& filenames);

  /**
    * @brief reads through all records in a stream (e.g. stdin) parsing them into Fastq
    * objects
    *
    * @param input a reference to the input stream (e.g. &std::cin)
    */
  explicit FastqReader(std::istream* const input);

  /**
    * @brief move constructor for the FastqReader class simply transfers all objects with the state
    * maintained.
    */
  FastqReader(FastqReader&&) = default;
  FastqReader& operator=(FastqReader&&) = default;

  /**
    * @brief a FastqReader cannot be copied safely, as it is iterating over a stream.
    */
  FastqReader(const FastqReader&) = delete;
  FastqReader& operator=(const FastqReader&) = delete;

  /**
    * @brief creates a FastqIterator pointing at the start of the input stream (needed by for-each
    * loop)
    *
    * @return a FastqIterator ready to start parsing the file
    */
  FastqIterator begin();

  /**
    * @brief creates a FastqIterator with a nullified input stream (needed by for-each loop)
    *
    * @return a FastqIterator that will match the end status of the iterator at the end of the stream
    */
  FastqIterator end();

private:
  std::shared_ptr<std::istream> m_input_stream; ///< a pointer to the input stream
};

}  // end of namespace

#endif // gamgee__fastq_reader__guard
