#ifndef __gamgee__variant_field_utils__
#define __gamgee__variant_field_utils__

#include <string>
#include <stdexcept>
#include <cstring> // for memcpy

#include <iostream>

namespace gamgee {
namespace utils {

/**
 * @brief an enumeration of the types in htslib for the format field values
 * @note these must match the order in the htslib defines in htslib/vcf.h
 */
enum class FormatFieldType {NIL = 0, INT8 = 1, INT16 = 2, INT32 = 3, FLOAT = 5, STRING = 7};


/**
  * @brief converts the value in an index from the byte array into VALUE_TYPE
  * 
  * The byte array's underlying data representation is record specific, meaning that even numbers (like Integer)
  * can be represented multiple ways across registers (some with uint8_t others with uint32_t...) dictated by
  * the the maximum value in the field.
  * 
  * This member function provides correct index location and appropriately creates a new value of 
  * VALUE_TYPE to return to the user. 
  *
  * @exception Will throw std::invalid_argument exception if trying to create a string out of a numeric format or 
  * vice-versa. All numeric type conversions are internally truncated or expanded accordingly.
  */
inline int32_t convert_data_to_integer(const uint8_t* data_ptr, const int index, const uint8_t num_bytes_per_value, const FormatFieldType& type) {
  const auto p = data_ptr + (index * num_bytes_per_value);
  auto v = 0;
  switch (type) {
    case FormatFieldType::NIL:
      return 0; // null returns zero -- seems reasonable to me
    case FormatFieldType::INT8:
    case FormatFieldType::INT16:
    case FormatFieldType::INT32:
    case FormatFieldType::FLOAT:
      memcpy(&v, p, num_bytes_per_value); 
      return int32_t(v); // truncate float if necessary
    case FormatFieldType::STRING:
      throw std::invalid_argument("user requested an integer value but underlying type is actually a string");
    default:
      return 0; // undefined type -- impossible to happen just so the compiler doesn't warn us 
  }
}

/**
 * @copydoc convert_data_to_integer(int, const uint8_t, const FormatFieldType&)
 */
inline float convert_data_to_float(const uint8_t* data_ptr, int index, const uint8_t num_bytes_per_value, const FormatFieldType& type) {
  const auto p = data_ptr + (index * num_bytes_per_value);
  switch (type) {
    case FormatFieldType::NIL:
      return 0.0; // null returns zero -- seems reasonable to me
    case FormatFieldType::INT8:
    case FormatFieldType::INT16:
    case FormatFieldType::INT32:
    case FormatFieldType::FLOAT:
      return *(reinterpret_cast<const float*>(p)); 
    case FormatFieldType::STRING:
      throw std::invalid_argument("user requested a float value but underlying type is actually a string");
    default:
      return 0; // undefined type -- impossible to happen just so the compiler doesn't warn us 
  }
}

/**
 * @copydoc convert_data_to_integer(int, const uint8_t, const FormatFieldType&)
 */
inline std::string convert_data_to_string(const uint8_t* data_ptr, int index, const uint8_t num_bytes_per_value, const FormatFieldType& type) {
  const auto p = data_ptr + (index * num_bytes_per_value);
  switch (type) {
    case FormatFieldType::NIL:
      return std::string{}; // null returns empty string -- seems reasonable to me
    case FormatFieldType::INT8:
      return std::to_string(*(reinterpret_cast<const int8_t*>(p))); // convert from integer to string if necessary
    case FormatFieldType::INT16:
      return std::to_string(*(reinterpret_cast<const int16_t*>(p))); // convert from integer to string if necessary
    case FormatFieldType::INT32:
      return std::to_string(*(reinterpret_cast<const int32_t*>(p))); // convert from integer to string if necessary
    case FormatFieldType::FLOAT:
      return std::to_string(*(reinterpret_cast<const float*>(p))); // convert from float to string if necessary
    case FormatFieldType::STRING:
      std::clog << "DEBUG: " << reinterpret_cast<const char *>(p) << std::endl;
      return std::string{reinterpret_cast<const char *>(p)}; // it is a simple string
    default:
      return std::string{}; // undefined type -- impossible to happen just so the compiler doesn't warn us 
  }
}

/**
 * @brief returns the number of bytes for a given FormatFieldType
 */
inline uint8_t size_for_type(const FormatFieldType& type) {
  switch (type) {
    case FormatFieldType::NIL:
    case FormatFieldType::INT8: 
    case FormatFieldType::INT16:
      return static_cast<uint8_t>(type);
    case FormatFieldType::INT32:
    case FormatFieldType::FLOAT:
      return 4;
    case FormatFieldType::STRING:
      return 1;
    default: 
      return 0; // 
  }
}

} // end namespace utils
} // end namespace gamgee

#endif
