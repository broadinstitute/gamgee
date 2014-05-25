#ifndef __gamgee__format_field_generic_value__
#define __gamgee__format_field_generic_value__

#include "variant_field_value_iterator.h"
#include "utils/hts_memory.h"
#include "utils/utils.h"

#include "htslib/vcf.h"

#include <memory>

namespace gamgee {

/**
 * @brief an enumeration of the types in htslib for the format field values
 * @note these must match the order in the htslib defines in htslib/vcf.h
 */
enum class FormatFieldType {NIL = 0, INT8 = 1, INT16 = 2, INT32 = 3, FLOAT = 5, STRING = 7};

/**
 * @brief A class template to hold the values a format field for a particular sample
 *
 * The role of this class is to perform the pointer manipulations behind the
 * scenes that permit the user to navigate the values of a field in a sample
 * without making any copies and benefiting from data locality (all the data is
 * stored contiguously in memory). 
 *
 * A typical use of the VariantFieldValue can be examplified by the
 * genotype qualitiy accessor in Variant: 
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * const auto all_pls = record.phred_likelihoods();
 * const auto my_pls = all_pls[2] // take the second sample
 * for_each(my_pls.begin(), my_pls.end(), [](const auto pl) { cout << pl << endl; });
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * VariantFieldValue objects can also be used in for loops like so: 
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (const auto pl : my_pls) 
 *   cout << pl << endl;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * While the VariantFieldValue objects are not really intended to be
 * created by the user, they are the returned by the FormatField iterator for
 * types that don't have a specialized object.  Utilizing them correctly can
 * really simplify your work by leveraging the power of the STL functions.
 *
 * @note all methods are inlined on purpose because they are so simple
 * 
 * @tparam TYPE the object type that holds the values for each sample. For
 * example for GQ it's a uint8_t, some types like GT and PL can have
 * specialized classes like Genotype and PhredLikelihoods. For all other types
 * it can be the GenericFormatFieldValue.
 */
template<class VALUE_TYPE>
class VariantFieldValue {
 public:

  /**
   * @brief creates a new VariantFieldValue poiinting to the shared byte array inside the variant object
   * @copydetails FormatField::FormatField(const std::shared_ptr<bcf1_t>&, bcf_fmt_t*)
   * @param body the the bcf1_t structure to hold a shared pointer to
   * @param format_ptr the format field pointer inside the body
   * @param data_ptr the location in the specific value inside the format_ptr byte array
   */
  VariantFieldValue(const std::shared_ptr<bcf1_t>& body, const bcf_fmt_t* const format_ptr, uint8_t* const data_ptr) :
    m_body {body},
    m_format_ptr {format_ptr},
    m_data_ptr {data_ptr},
    m_num_bytes {size_for_type(format_ptr->type)}
  {}

  /**
   * @brief copying of the VariantFieldValue object is not allowed.
   * @copydetails FormatField::FormatField(const FormatField&)
   */
  VariantFieldValue(const VariantFieldValue& other) = delete;
  
  /**
   * @brief safely moves the data from one VariantFieldValue to a new one without making any copies
   * @param other another VariantFieldValue object
   */
  VariantFieldValue(VariantFieldValue&& other) :
    m_body {std::move(other.m_body)},
    m_format_ptr {other.m_format_ptr},
    m_data_ptr {other.m_data_ptr},
    m_num_bytes {other.m_num_bytes}
  {}

  /**
   * @copydoc FormatField::VariantFieldValue(const VariantFieldValue&)
   */
  VariantFieldValue& operator=(const VariantFieldValue& other) = delete;

  /**
   * @brief safely moves the data from one FormatField to the other without making any copies
   * @param other another VariantFieldValue object
   */
  VariantFieldValue& operator=(VariantFieldValue&& other) {
    if (this != &other)
      return *this;
    m_body = std::move(other.m_body);
    m_format_ptr = other.m_format_ptr;
    m_data_ptr = other.m_data_ptr;
    m_num_bytes = other.m_num_bytes;
    return *this;
  }

  /**
   * @brief random access to the value of a given sample for reading or writing
   * @param index must be between 0 and the number of values per sample for this record 
   * @note implementation guarantees this operation to be O(1)
   * @exception std::out_of_range if sample is out of range
   * @return the value in that index
   */
  VALUE_TYPE operator[](const uint32_t index) const {
    utils::check_boundaries(index, m_format_ptr->n - 1);
    return convert_from_byte_array(index); 
  }

  VariantFieldValueIterator<VALUE_TYPE> begin() const {
    return VariantFieldValueIterator<VALUE_TYPE>{m_body, m_data_ptr};
  }

  VariantFieldValueIterator<VALUE_TYPE> end() const {
    return VariantFieldValueIterator<VALUE_TYPE>{m_body, m_data_ptr + m_format_ptr->size};
  }

 private:
  const std::shared_ptr<bcf1_t> m_body;
  const bcf_fmt_t* const m_format_ptr;
  uint8_t* const m_data_ptr;
  const uint8_t m_num_bytes;

  VALUE_TYPE convert_from_byte_array(int index) const;

  static uint8_t size_for_type(int type) {
    switch (type) {
      case 0:
      case 1: 
      case 2:
        return type;
      case 3:
      case 5:
        return 4;
      case 7:
        return 1;
      default: 
        return 0; // 
    }
  }
};

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

template<> inline
int32_t VariantFieldValue<int32_t>::convert_from_byte_array(int index) const {
  const auto p = m_data_ptr + (index * m_num_bytes);
  auto v = 0;
  switch (static_cast<FormatFieldType>(m_format_ptr->type)) {
    case FormatFieldType::NIL:
      return 0; // null returns zero -- seems reasonable to me
    case FormatFieldType::INT8:
    case FormatFieldType::INT16:
    case FormatFieldType::INT32:
    case FormatFieldType::FLOAT:
      memcpy(&v, p, m_num_bytes); 
      return int32_t(v); // truncate float if necessary
    case FormatFieldType::STRING:
      throw std::invalid_argument("user requested an integer value but underlying type is actually a string");
    default:
      return 0; // undefined type -- impossible to happen just so the compiler doesn't warn us 
  }
}

/**
 * @copydoc VariantFieldValue<int32_t>::convert_from_byte_array(int)
 */
template<> inline
float VariantFieldValue<float>::convert_from_byte_array(int index) const {
  const auto p = m_data_ptr + (index * m_num_bytes);
  auto v = 0.0;
  switch (static_cast<FormatFieldType>(m_format_ptr->type)) {
    case FormatFieldType::NIL:
      return 0.0; // null returns zero -- seems reasonable to me
    case FormatFieldType::INT8:
    case FormatFieldType::INT16:
    case FormatFieldType::INT32:
    case FormatFieldType::FLOAT:
      memcpy(&v, p, m_num_bytes); 
      return float(v); // convert from integer to float if necessary
    case FormatFieldType::STRING:
      throw std::invalid_argument("user requested a float value but underlying type is actually a string");
    default:
      return 0; // undefined type -- impossible to happen just so the compiler doesn't warn us 
  }
}

/**
 * @copydoc VariantFieldValue<int32_t>::convert_from_byte_array(int)
 */
template<> inline
std::string VariantFieldValue<std::string>::convert_from_byte_array(int index) const {
  const auto p = m_data_ptr + (index * m_num_bytes);
  auto v = 0;
  switch (static_cast<FormatFieldType>(m_format_ptr->type)) {
    case FormatFieldType::NIL:
      return std::string{}; // null returns empty string -- seems reasonable to me
    case FormatFieldType::INT8:
    case FormatFieldType::INT16:
    case FormatFieldType::INT32:
    case FormatFieldType::FLOAT:
      memcpy(&v, p, m_num_bytes); 
      return std::to_string(v);
    case FormatFieldType::STRING:
      return std::string{reinterpret_cast<char *>(p)};
    default:
      return std::string{}; // undefined type -- impossible to happen just so the compiler doesn't warn us 
  }
}


}

#endif
