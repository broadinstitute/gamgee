#ifndef gamgee__fastq_reader__guard
#define gamgee__fastq_reader__guard

#include "fastq_iterator.h"

#include <string>
#include <iostream>
#include <fstream>

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
    * @brief reades through all records in a file (fasta or fastq) parsing them into Fastq
    * objects
    *
    * @param filename the name of the fasta/fastq file
    */
  FastqReader(const std::string& filename);

  /**
    * @brief reades through all records in a stream (e.g. stdin) parsing them into Fastq
    * objects
    *
    * @param input a reference to the input stream (e.g. &std::cin)
    */
  FastqReader(std::istream* const input);

  /**
    * @brief closes the file stream if there is one (in case we are reading a fasta/fastq file)
    */
  ~FastqReader();

  /**
    * @brief move constructor for the FastqReader class simply transfers all objects with the state
    * maintained.
    */
  FastqReader(FastqReader&&);

  /**
    * @brief a FastqReader cannot be copied safely, as it is iterating over a stream.
    */
  FastqReader(const FastqReader&) = delete;

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
  std::istream* m_input_stream;  ///< a pointer to the input stream
  std::ifstream m_file_stream;   ///< a file stream (in case we are reading from a file)
};

}  // end of namespace

#endif // gamgee__fastq_reader__guard
