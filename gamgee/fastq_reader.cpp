#include "fastq_reader.h"
#include "fastq_iterator.h"

#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <limits>

using namespace std;

namespace gamgee {

FastqReader::FastqReader(const std::string& filename) {
  if (!filename.empty()) {
    m_file_stream.open(filename);
    m_input_stream = &m_file_stream;
  }
}

FastqReader::~FastqReader() {
  m_file_stream.close();
}

FastqReader::FastqReader(std::istream* const input) : 
  m_input_stream{input} 
{}

FastqIterator FastqReader::begin() {
  return FastqIterator{m_input_stream};
}

FastqIterator FastqReader::end() {
  return FastqIterator{};
}


}  // end of namespace
