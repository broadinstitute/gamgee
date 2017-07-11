#ifndef gamgee__variant_field_utils__guard
#define gamgee__variant_field_utils__guard

#include "htslib/vcf.h"

#include <string>
#include <iostream>

namespace gamgee {
namespace utils {

/**
 * @brief an enumeration of the types in htslib for the format field values
 * @note these must match the order in the htslib defines in htslib/vcf.h
 */
enum class VariantFieldType {NIL = 0, INT8 = 1, INT16 = 2, INT32 = 3, FLOAT = 5, STRING = 7};


/**
  * @brief converts the value in an index from the byte array into int32_t
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
int32_t convert_data_to_integer(const uint8_t* data_ptr, const int index, const uint32_t num_bytes_per_value, const VariantFieldType& type);

/**
 * @brief converts the value in an index from the byte array into float
 * @copydetails convert_data_to_integer(int, const uint8_t, const VariantFieldType&)
 */
float convert_data_to_float(const uint8_t* data_ptr, const int index, const uint32_t num_bytes_per_value, const VariantFieldType& type);

/**
 * @brief converts the value in an index from the byte array into string
 * @copydetails convert_data_to_integer(int, const uint8_t, const VariantFieldType&)
 */
std::string convert_data_to_string(const uint8_t* data_ptr, const int index, const uint32_t num_bytes_per_value, const VariantFieldType& type);

/**
 * @brief returns the number of bytes for a given VariantFieldType
 */
uint32_t size_for_type(const VariantFieldType& type, const bcf_fmt_t* const format_ptr);

/**
 * @brief returns the number of bytes for a given VariantFieldType
 */
uint32_t size_for_type(const VariantFieldType& type, const bcf_info_t* const info_ptr);

/**
 * @brief - check if type is of type string
 */
inline bool is_string_type(const int32_t& type) {
  return (static_cast<utils::VariantFieldType>(type) == VariantFieldType::STRING); 
}

} // end namespace utils
} // end namespace gamgee

#endif // gamgee__variant_field_utils__guard
