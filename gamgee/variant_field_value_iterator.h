#ifndef __gamgee__variant_field_value_iterator__
#define __gamgee__variant_field_value_iterator__

#include "utils/format_field_type.h"

#include "htslib/vcf.h"

#include <iterator>
#include <memory>

namespace gamgee {

/** 
 * @brief iterator for FormatFieldGenericValue objects. 
 * 
 * This iterator will walk through all the values for a sample for a given
 * FormatFieldGenericValue object. For example if you want to iterate over all
 * the GQ values of a Variant record you would do so through this iterator.
 *
 * @note implements a random access iterator which gives you full performance
 * on STL algorithms that use iterators (mostly every one)
 * 
 * @note this iterator never makes any copies of the underlying memory. It is
 * implemented with pointer acrobatics behind the scenes for maximum
 * performance while maintaining a friendly interface.
 *
 * @warning modifying any elements via this iterator **will modify** the values
 * in the Variant record. 
 */
template<class VALUE_TYPE>
class VariantFieldValueIterator : public std::iterator<std::random_access_iterator_tag, VALUE_TYPE> {
 public:

  /**
   * @brief simple constructor
   * @param body a pointer to the bcf1_t data structure to be held as a shared pointer
   * @param data_ptr the byte array containing the value(s) for this sample
   * @note this constructor is probably only used by VariantFieldValue::begin() and
   * VariantFieldValue::end()
   */
  VariantFieldValueIterator(const std::shared_ptr<bcf1_t>& body, uint8_t* data_ptr, const uint8_t num_bytes, const utils::FormatFieldType& type) :
    m_body {body},
    m_current_data_ptr {data_ptr},
    m_original_data_ptr {data_ptr},
    m_num_bytes {num_bytes},
    m_type {type}
  {}
      
  /**
   * @copydoc VariantFieldIterator::VariantFieldIterator(const VariantFieldIterator&)
   */
  VariantFieldValueIterator(const VariantFieldValueIterator& other) :
    m_body {other.m_body},
    m_current_data_ptr {other.m_current_data_ptr},
    m_original_data_ptr {other.m_original_data_ptr},
    m_num_bytes {other.m_num_bytes},
    m_type {other.m_type}
  {}

  /**
   * @copydoc VariantFieldIterator::VariantFieldIterator(const VariantFieldIterator&)
   */
  VariantFieldValueIterator(VariantFieldValueIterator&& other) noexcept :
    m_body {std::move(other.m_body)},
    m_current_data_ptr {other.m_current_data_ptr},
    m_original_data_ptr {other.m_original_data_ptr},
    m_num_bytes {other.m_num_bytes},
    m_type {other.m_type}
  {}

  /**
   * @copydoc VariantFieldIterator::VariantFieldIterator(const VariantFieldIterator&)
   */
  VariantFieldValueIterator& operator=(const VariantFieldValueIterator& other) {
    if (&this == other)
      return *this;
    m_body = std::move(other.m_body);
    m_current_data_ptr = other.m_current_data_ptr;
    m_original_data_ptr = other.m_original_data_ptr;
    m_num_bytes = other.m_num_bytes;
    m_type = other.m_type;
    return *this;
  }


  /**
   * @copydoc VariantFieldIterator::VariantFieldIterator(VariantFieldIterator&&)
   */
  VariantFieldValueIterator& operator=(VariantFieldValueIterator&& other) noexcept {
    if (&this == other)
      return *this;
    m_body = std::move(other.m_body);
    m_current_data_ptr = other.m_current_data_ptr;
    m_original_data_ptr = other.m_original_data_ptr;
    m_num_bytes = other.m_num_bytes;
    m_type = other.m_type;
    return *this;
  }

  /**
   * @copydoc VariantFieldIterator::operator+=(const int)
   */
  VariantFieldValueIterator& operator+=(const int n) {
    m_current_data_ptr += n * m_num_bytes;
    return *this;
  }

  /**
   * @copydoc VariantFieldIterator::operator+=(int)
   */
  VariantFieldValueIterator& operator-=(const int n) {
    m_current_data_ptr -= n * m_num_bytes;
    return *this;
  }

  /**
   * @brief two iterators are equal if they are in exactly the same state (pointing at the same location in memory
   */
  bool operator==(const VariantFieldValueIterator& other) {
    return m_body == other.m_body && m_current_data_ptr == other.m_current_data_ptr;
  }

  /**
   * @brief the oposite check of VariantFieldValueIterator::operator==()
   */
  bool operator!=(const VariantFieldValueIterator& other) {
    return m_body != other.m_body || m_current_data_ptr != other.m_current_data_ptr;
  }

  /**
   * @brief an operator is greater/less than another iterator if it is pointing to a previous element (sample) in the FormatField 
   * object. The order is determined by the Variant record.
   */
  bool operator<(const VariantFieldValueIterator& other) {
    return m_body == other.m_body && m_current_data_ptr < other.m_current_data_ptr;
  }

  /**
   * @copydoc VariantFieldValueIterator::operator<()
   */
  bool operator>(const VariantFieldValueIterator& other) {
    return m_body == other.m_body && m_current_data_ptr > other.m_current_data_ptr;
  }

  /**
   * @copydoc VariantFieldValueIterator::operator<()
   */
  bool operator<=(const VariantFieldValueIterator& other) {
    return m_body == other.m_body && m_current_data_ptr <= other.m_current_data_ptr;
  }

  /**
   * @copydoc VariantFieldValueIterator::operator<()
   */
  bool operator>=(const VariantFieldValueIterator& other) {
    return m_body == other.m_body && m_current_data_ptr >= other.m_current_data_ptr;
  }

  /**
   * @brief direct access to the value of the current sample
   * @return the value in its native type
   */
  VALUE_TYPE operator*() const noexcept {
    return convert_from_byte_array(m_current_data_ptr, 0);
  }

  /**
   * @brief advances to the next sample
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result of operator* with end(). This is the STL way.
   * @return the next value in it's native type
   */
  VALUE_TYPE operator++() noexcept {
    m_current_data_ptr += m_num_bytes;
    return convert_from_byte_array(m_current_data_ptr, 0);
  }

  /**
   * @brief advances to the previous sample
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result of operator* with end(). This is the STL way.
   * @return the previous value in it's native type
   */
  VALUE_TYPE operator--() {
    m_current_data_ptr -= m_num_bytes;
    return convert_from_byte_array(m_current_data_ptr, 0);
  }

  /**
   * @brief random access to the value of a given index for reading or writing
   * @param index must be between 0 and the number of indices for this record but no boundary check is done in this implementation
   * @note implementation guarantees this operation to be O(1)
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result of operator* with end(). This is the STL way.
   * @return the value in it's native type
   */
  VALUE_TYPE operator[](const uint32_t index) const {
    return convert_from_byte_array(m_original_data_ptr, index);
  }

 private:
  const std::shared_ptr<bcf1_t> m_body;
  const uint8_t* m_current_data_ptr;
  const uint8_t* const m_original_data_ptr;
  const uint8_t m_num_bytes;
  const utils::FormatFieldType m_type;

  VALUE_TYPE convert_from_byte_array(const uint8_t* data_ptr, int index) const;
};

/**
 * @brief specialization of the conversion template for int32_t
 */
template<> inline
int32_t VariantFieldValueIterator<int32_t>::convert_from_byte_array(const uint8_t* data_ptr, int index) const {
  return utils::convert_data_to_integer(data_ptr, index, m_num_bytes, m_type);
}

/**
 * @brief specialization of the conversion template for floats
 */
template<> inline
float VariantFieldValueIterator<float>::convert_from_byte_array(const uint8_t* data_ptr, int index) const {
  return utils::convert_data_to_float(data_ptr, index, m_num_bytes, m_type);
}

/**
 * @brief specialization of the conversion template for strings
 */
template<> inline
std::string VariantFieldValueIterator<std::string>::convert_from_byte_array(const uint8_t* data_ptr, int index) const {
  return utils::convert_data_to_string(data_ptr, index, m_num_bytes, m_type);
}

};

#endif
