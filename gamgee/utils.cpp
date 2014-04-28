#include "utils.h" 

#include <string>
#include <algorithm>

namespace gamgee {
namespace utils {

char complement_base (const char base) {
  switch (base) {
    case 'A':
      return 'T';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    case 'T':
      return 'A';
    case 'a':
      return 't';
    case 'c':
      return 'g';
    case 'g':
      return 'c';
    case 't':
      return 'a';
    default:
      return base;
  }
}

std::string complement(std::string& seq) {
  std::transform(seq.begin(), seq.end(), seq.begin(), complement_base);
  return seq;
}

std::string complement(const std::string& seq) {
  auto rev = std::string{};
  rev.reserve(seq.length());
  std::transform(seq.cbegin(), seq.cend(), std::back_inserter(rev), complement_base);
  return rev;
}

std::string reverse_complement(const std::string& seq) {
  auto rev = std::string{};
  rev.reserve(seq.length());
  std::transform(seq.crbegin(), seq.crend(), std::back_inserter(rev), complement_base);
  return rev;
}

}
}
