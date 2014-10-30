#include "fastq_reader.h"
#include "fastq_iterator.h"

#include "exceptions.h"
#include "utils/file_utils.h"

#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <limits>

using namespace std;

namespace gamgee {

FastqReader::FastqReader(const std::string& filename) :
  m_input_stream {}
{
  if (!filename.empty()) {
    init_reader(filename);
  }
}

FastqReader::FastqReader(const std::vector<std::string>& filenames) :
  m_input_stream {}
{
  if (filenames.size() > 1)
    throw SingleInputException{"filenames", filenames.size()};
  if (!filenames.empty()) {
    init_reader(filenames.front());
  }
}

FastqReader::FastqReader(std::istream* const input) : 
  m_input_stream{shared_ptr<std::istream>(input)}
{}

FastqIterator FastqReader::begin() {
  return FastqIterator{m_input_stream};
}

FastqIterator FastqReader::end() {
  return FastqIterator{};
}

void FastqReader::init_reader(const std::string& filename) {
  m_input_stream = utils::make_shared_ifstream(filename);
  if ( m_input_stream->fail() ) {
    throw FileOpenException{filename};
  }
}


}  // end of namespace
