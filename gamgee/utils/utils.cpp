#include "utils.h" 

#include <string>
#include <algorithm>
#include <memory.h>
#include <vector>

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

std::string complement(std::string& sequence) {
  std::transform(sequence.begin(), sequence.end(), sequence.begin(), complement_base);
  return sequence;
}

std::string complement(const std::string& sequence) {
  auto rev = std::string{};
  rev.reserve(sequence.length());
  std::transform(sequence.cbegin(), sequence.cend(), std::back_inserter(rev), complement_base);
  return rev;
}

std::string reverse_complement(const std::string& sequence) {
  auto rev = std::string{};
  rev.reserve(sequence.length());
  std::transform(sequence.crbegin(), sequence.crend(), std::back_inserter(rev), complement_base);
  return rev;
}

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) { 
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...)); 
}

std::vector<std::string> hts_string_array_to_vector(const char * const * const string_array, const uint32_t array_size) {
  auto result = std::vector<std::string>{};
  result.reserve(array_size);
  for (auto i = 0u; i != array_size; ++i) 
    result.emplace_back(string_array[i]);
  return result;
}

void check_boundaries(const int index, const int max_index, const int min_index) {
  if (index < min_index || index > max_index) {
    std::stringstream error_message {};  ///< @todo update this to auto when gcc-4.9 is available on travis-ci
    error_message << "The index requested is out of range: " << index << " the maximum index is " << max_index << " and the minimum is " << min_index << std::endl;
    throw std::out_of_range(error_message.str());
  }
}


}
}
