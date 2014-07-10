#ifndef __gamgee_utils__
#define __gamgee_utils__

#include <string>
#include <memory>
#include <vector>
#include <sstream>

namespace gamgee {

/** 
 * @brief utility functions for the gamgee library
 *
 * This namespace includes typical functors and templates for DNA sequence manipulation as
 * well as low level routines for htslib memory management
 */
namespace utils {

/**
 * @brief calculates the reverse complement of a sequence
 * @param sequence the sequence to get the reverse complement from
 * @return a new string with the reverse complement of the sequence
 */
std::string reverse_complement(const std::string& sequence);

/**
 * @brief calculates the complement of a sequence
 * @param sequence the sequence to get the complement from
 * @return a new string with the complement of the sequence
 */
std::string complement(const std::string& sequence); 

/**
 * @brief calculates the complement of a sequence in-place
 * @param sequence the sequence to turn into the complement
 * @return the complement of the sequence calculated in-place (modifies the sequence)
 */
std::string complement(std::string& sequence);               

/**
 * @brief calculates the complement of a base
 * @param base the base to get the complement from
 * @return a new char with the complement of base
 */
char complement (const char base);

/**
 * @brief herb sutter's implementation of make unique
 * @note this is going to become standard in C++14, we can drop this as soon as gcc 4.9 and clang 3.5 are released
 */
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args);

/**
 * @brief converts an array of c-strings into a vector<string>
 * 
 * Useful in many hts structs where a list of names is stored as a char** and we want to 
 * manipulate it in a vector<string> to maintain data contiguity and improve usability
 */
std::vector<std::string> hts_string_array_to_vector(const char * const * const string_array, const uint32_t array_size);

/**
 * @brief checks that an index is within min_index to max_index
 * @param index the index to check
 * @param max_index the maximum valid index
 * @param min_index the minimum valid index (typically 0 for array boundaries, but can be arbitrary for pointer checks)
 * @exception throws an out_of_bounds exception if index is out of limits
 */
inline void check_boundaries(const int index, const int max_index, const int min_index = 0) {
  if (index < min_index || index > max_index) {
    std::stringstream error_message {};  ///< @todo update this to auto when gcc-4.9 is available on travis-ci
    error_message << "The index requested is out of range: " << index << " the maximum index is " << max_index << " and the minimum is " << min_index << std::endl;
    throw std::out_of_range(error_message.str());
  }
}


/**
 * @brief checks that an index is greater than max_index
 * @param index the index to check
 * @param max_index the maximum valid index
 * @exception throws an out_of_bounds exception if index is out of limits
 */
inline void check_max_boundary(const uint32_t index, const uint32_t max_index) {
  if (index > max_index) {
    std::stringstream error_message {};  ///< @todo update this to auto when gcc-4.9 is available on travis-ci
    error_message << "The index requested is out of range: " << index << " the maximum index is " << max_index << std::endl;
    throw std::out_of_range(error_message.str());
  }
}


} // end utils namespace
} // end gamgee namespace

#endif
