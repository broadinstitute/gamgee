#ifndef gamgee__utils__guard
#define gamgee__utils__guard

#include <string>
#include <memory>
#include <vector>
#include <sstream>
#include "htslib/vcf.h"

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
 * @brief converts an array of c-strings into a vector<string>
 * 
 * Useful in many hts structs where a list of names is stored as a char** and we want to 
 * manipulate it in a vector<string> to maintain data contiguity and improve usability
 */
std::vector<std::string> hts_string_array_to_vector(const char * const * const string_array, const uint32_t array_size);

/**
 * @brief checks that an index is greater than or equal to size
 * @param index the index between 0 and size to check
 * @param size one past the maximum valid index
 * @param prefix_msg additional string to prefix error message
 * @exception throws an out_of_bounds exception if index is out of limits
 */
inline void check_max_boundary(const uint32_t index, const uint32_t size, const std::string& prefix_msg) {
  if (index >= size) {
    std::stringstream error_message {};  ///< @todo update this to auto when gcc-4.9 is available on travis-ci
    error_message << prefix_msg << "Index:  " << index << " must be less than " << size << std::endl;
    throw std::out_of_range(error_message.str());
  }
}

/**
 * @brief checks that an index is greater than or equal to size
 * @param index the index between 0 and size to check
 * @param size one past the maximum valid index
 * @exception throws an out_of_bounds exception if index is out of limits
 */
inline void check_max_boundary(const uint32_t index, const uint32_t size) {
  if (index >= size) {
    std::stringstream error_message {};  ///< @todo update this to auto when gcc-4.9 is available on travis-ci
    error_message << "Index:  " << index << " must be less than " << size << std::endl;
    throw std::out_of_range(error_message.str());
  }
}
/**
 * @brief - Check whether two values from VCF fields of primitive types (for which the == operator is defined) * are equal
 * The function template is specialized for float fields
 */
template<class TYPE>
inline bool bcf_check_equal_element(const TYPE& x, const TYPE& y) {
  return (x == y);
}
/**
 * @brief: Check whether two float values from VCF fields  are equal
 * @note: since bcf_float_missing and bcf_float_vector_end are represented as NaNs and (nan == nan) returns false,
 * a special check needs to be done for checking equality of float values in VCFs
 */
template<>
inline bool bcf_check_equal_element<float>(const float& x, const float& y) {
  return ((x == y) || (bcf_float_is_missing(x) && bcf_float_is_missing(y))
	|| (bcf_float_is_vector_end(x) && bcf_float_is_vector_end(y)));
}

/*
 * @brief: returns true if value is the vector end for current TYPE
 * @note: for non numeric types returns false, as vector end is undefined
 */
template<class TYPE> inline
bool bcf_is_vector_end_value(const TYPE& value) {
  return false;
}

/*
 * @brief specialization of the bcf_is_vector_end_value function for int32_t
 */
template<> inline
bool bcf_is_vector_end_value<int32_t>(const int32_t& value) {
  return (value == bcf_int32_vector_end);
}

/*
 * @brief specialization of the bcf_is_vector_end_value function for float
 */
template<> inline
bool bcf_is_vector_end_value<float>(const float& value) {
  return bcf_float_is_vector_end(value);
}



/**
 * @brief advances current ptr to end of the vector if the current element is bcf_*_vector_end
 */
template<class ITER> inline
const uint8_t* cache_and_advance_to_end_if_necessary(const uint8_t* current_ptr, const uint8_t* end_ptr, ITER& it) {
  if(end_ptr != nullptr && current_ptr < end_ptr 
      && utils::bcf_is_vector_end_value(it.read_and_cache_current_pointee()))
    return end_ptr;
  return current_ptr;
}

/**
 * A special version of starts_with for specific tag search.
 *
 * Shared utility function for ReadGroup and Program classes to parse SamHeader
 * raw input. Useful to check if the result of std::string::find caught a token
 * of interest.
 */
inline bool starts_with(const std::string& token, const char tag[2]) {
    return token[0] == tag[0] && token[1] == tag[1];
}

/**
 * @brief Store the input integer value in little-endian byte order
 */
template<class TYPE>
inline uint8_t* to_little_endian(const TYPE& value, uint8_t *buffer) {
  for (std::size_t i = 0 ; i < sizeof(TYPE) ; ++i) {
      *buffer++ = (uint8_t) ((value >> (i * 8)) & 0xff);
   }
   return buffer;
}




} // end utils namespace
} // end gamgee namespace

#endif // gamgee__utils__guard
