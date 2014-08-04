#include <string>
#include <iostream>
#include <limits>

#include "fastq_iterator.h"

using namespace std;

namespace gamgee {

FastqIterator::FastqIterator() {
  m_input_stream = nullptr;
}

FastqIterator::FastqIterator(std::istream* in) : 
  m_input_stream {in} 
{
  m_is_fastq = false;
  if (m_input_stream->get() == '@') {
    m_is_fastq  = true;
    m_bor_delim = '@';
    m_eos_delim = '+';
  } else {
    m_is_fastq  = false;
    m_bor_delim = '>';
    m_eos_delim = '>';
  }
  m_element = fetch_next_element();
}

Fastq& FastqIterator::operator*() {
  return m_element;
}

Fastq& FastqIterator::operator++() {
  m_element = fetch_next_element();
  return m_element;
}

bool FastqIterator::operator!=(const FastqIterator& rhs) {
  return m_input_stream != rhs.m_input_stream;
}

const string FastqIterator::parse_multiline() {
  auto s = string{};
  *m_input_stream >> s; 
  skip_new_lines();
  return s;
}

string FastqIterator::parse_comment() {
  auto comment = string{};
  while (m_input_stream->peek() == ' ')
    m_input_stream->ignore();
  getline(*m_input_stream, comment);
  return comment;
}

string FastqIterator::parse_seq() {
  auto seq = string{};
  while(m_input_stream->good() && m_input_stream->peek() != m_eos_delim) {
    seq += parse_multiline();
  }
  return seq;
}

string FastqIterator::parse_quals(uint32_t seq_length) {
  auto qual = string{};
  if (m_is_fastq) {
    m_input_stream->ignore(numeric_limits<streamsize>::max(), '\n');  // skip the line with the second record name
    while (m_input_stream->good() && qual.length() < seq_length ) {
      qual += parse_multiline();
    }
  } 
  m_input_stream->ignore(numeric_limits<streamsize>::max(), m_bor_delim); //skip everything (even if fastq is malformed and has extra quals) until the next record. It's essential that we remove the '@' here.
  return qual;
}

Fastq FastqIterator::fetch_next_element() {
  auto name = string{};
  if (!(*m_input_stream >> name)) { // parse name (first string) and abort if we reached the end of the file
    m_input_stream = nullptr;
    return Fastq{};
  }
  const auto comment = parse_comment();
  const auto seq     = parse_seq();
  const auto quals   = parse_quals(seq.length());
  return Fastq {move(name), move(comment), move(seq), move(quals)};
}

void FastqIterator::skip_new_lines() {
  while (m_input_stream->peek() == '\n') 
    m_input_stream->ignore(numeric_limits<streamsize>::max(), '\n');
}

}

