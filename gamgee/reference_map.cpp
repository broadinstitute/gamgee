#include "reference_map.h"

#include "fastq_reader.h"
#include "interval.h"
#include "utils/utils.h"

#include <unordered_map>
#include <string>

using namespace gamgee;
using namespace std;

namespace gamgee {

ReferenceMap::ReferenceMap(const std::string& filename) {
  FastqReader reader{filename};
  read_fastq(reader);
}

void ReferenceMap::read_fastq(FastqReader& reader) {
  for (auto& fq : reader)  
    this->insert({fq.name(), fq.sequence()});
}

string ReferenceMap::get_sequence(const Interval& interval, const bool reverse_strand) const {
  const auto& seq = this->at(interval.chr());
  auto result = seq.substr(interval.start()-1, interval.size());
  return reverse_strand ? utils::complement(result) : result;
}

}  // namespace gamgee
