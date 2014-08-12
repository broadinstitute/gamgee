#include "fastq_reader.h"
#include "fastq_iterator.h"

#include "utils/utils.h"

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

FastqReader::FastqReader(const std::vector<std::string>& filenames) {
  if (filenames.size() > 1)
    throw utils::SingleInputException{"filenames", filenames.size()};
  if (!filenames.empty()) {
    m_file_stream.open(filenames.front());
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
