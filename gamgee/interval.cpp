#include "interval.h"

#include<boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <sstream>
#include <fstream>

using namespace std;
using namespace boost;

namespace gamgee {

const auto PICARD_HEADER_TAG = '@';
const auto sep = char_separator<char>{" \t:-"};

void skip_picard_header(istream& infile) {
  auto line = string{};
  while(infile.peek() == PICARD_HEADER_TAG)
    getline(infile, line);
}

Interval parse_interval_record(const string& line) {
  auto tokens = tokenizer<char_separator<char>>{line, sep};
  auto chr = string{};
  auto start = 0u;
  auto stop = -1;
  auto it = tokens.begin();
  chr = *it;
  if( ++it !=  tokens.end()) 
    start = lexical_cast<uint32_t>(*it);
  if( ++it != tokens.end()) 
    stop = lexical_cast<int32_t>(*it); // has to be int so we can check for failures
  if (stop == -1) 
    stop = start;  // interval may not have a stop (meaning one base interval). This makes sure stop is unsigned
  return Interval {chr, start, static_cast<uint32_t>(stop)};  // stop is guaranteed to be unsigned by this point
}

vector<Interval> read_intervals(const string& intervals_file) {
  ifstream infile {intervals_file};
  return read_intervals(infile);
}

vector<Interval> read_intervals(istream& input) {
  auto result = vector<Interval>{};
  skip_picard_header(input);
  auto line = string{};
  while(getline(input, line)) 
    result.emplace_back(parse_interval_record(line));
  return result;
}

vector<Interval> Interval::tsca_tiling(const uint32_t spacing, const uint32_t insert_size, const uint32_t flanking) {
  auto result = vector<Interval>{};
  const auto tile_length = insert_size + spacing;
  const auto start = max(1u, m_start - flanking);
  for (auto top_start = start; top_start + insert_size <= m_stop + flanking; top_start += tile_length) {
    const auto bottom_start = top_start + (insert_size+spacing)/2;  // start the bottom bait such that it is centered in the spacing between the top two baits
    result.emplace_back(m_chr, top_start, top_start+insert_size-1);
    result.emplace_back(m_chr, bottom_start , bottom_start+insert_size-1);
  }
  return result; 
}

vector<Interval> Interval::tile_left(const uint32_t tile_size, const uint32_t spacing) const {
  auto result = vector<Interval>{};
  const auto l = tile_size + spacing; // length of the jump to the next tile
  for(auto i = m_start; i+tile_size-1 <= m_stop; i += l) 
    result.emplace_back(m_chr, i, i+tile_size-1);
  return result;
}

vector<Interval> Interval::tile_right(const uint32_t tile_size, const uint32_t spacing) const {
  auto result = vector<Interval>{};
  const auto l = tile_size + spacing; // length of the jump to the next tile
  for(auto i = m_stop; i-tile_size+1 >= m_start; i -= l) 
    result.emplace_back(m_chr, i-tile_size+1, i);
  return result;
}

bool Interval::operator==(const Interval& rhs) const {
  return (this == &rhs) || ((m_chr == rhs.m_chr) && (m_start == rhs.m_start) && (m_stop == rhs.m_stop));
}

}  // end of namespace

ostream& operator<< (std::ostream& os, const gamgee::Interval& i) {
  using IT = gamgee::Interval::IntervalType;
  switch (i.output_type()) {
    case IT::PICARD:
      return os << i.chr() << "\t" << i.start() << "\t" << i.stop() << "\t+" ;
    case IT::BED:
      return os << i.chr() << "\t" << i.start() << "\t" << i.stop();
    default:  // GATK
      return os << i.chr() << ":" << i.start() << "-" << i.stop();
  }
}
