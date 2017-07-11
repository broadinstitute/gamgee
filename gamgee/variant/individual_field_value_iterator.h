#ifndef gamgee__individual_field_value_iterator__guard
#define gamgee__individual_field_value_iterator__guard

#include "../utils/variant_field_type.h"

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
class IndividualFieldValueIterator : public std::iterator<std::random_access_iterator_tag, VALUE_TYPE> {
 public:

   /**
   * @brief Constructor with bcf1_t structure and start and end pointers of the array/vector 
   * @param body a pointer to the bcf1_t data structure to be held as a shared pointer
   * @param data_ptr the byte array containing the value(s) for this sample
   * @param end_ptr pointer to the end of the byte array
   * @param num_bytes number of bytes in each value
   * @param type the encoding of the value 
   * @note this constructor is probably only used by IndividualFieldValue::begin() and
   * IndividualFieldValue::end()
   */
  explicit IndividualFieldValueIterator(const std::shared_ptr<bcf1_t>& body, uint8_t* data_ptr, uint8_t* end_ptr, const uint32_t num_bytes, const utils::VariantFieldType& type) :
    m_body {body},
    m_current_data_ptr {data_ptr},
    m_original_data_ptr {data_ptr},
    m_end_data_ptr {end_ptr},
    m_num_bytes {num_bytes},
    m_type {type},
    m_is_current_pointee_cached {false}
  {}

  /**
   * @brief Constructor with bcf1_t structure and start pointer of the array/vector 
   * @param body a pointer to the bcf1_t data structure to be held as a shared pointer
   * @param data_ptr the byte array containing the value(s) for this sample
   * @param num_bytes number of bytes in each value
   * @param type the encoding of the value 
   * @note this constructor is probably only used by IndividualFieldValue::begin() and
   * IndividualFieldValue::end()
   */
  explicit IndividualFieldValueIterator(const std::shared_ptr<bcf1_t>& body, uint8_t* data_ptr, const uint8_t num_bytes, const utils::VariantFieldType& type): IndividualFieldValueIterator(body, data_ptr, nullptr, num_bytes, type)
  {}
  
  /**
   * @copydoc IndividualFieldIterator::IndividualFieldIterator(const IndividualFieldIterator&)
   */
  IndividualFieldValueIterator(const IndividualFieldValueIterator& other) :
    m_body {other.m_body},
    m_current_data_ptr {other.m_current_data_ptr},
    m_original_data_ptr {other.m_original_data_ptr},
    m_end_data_ptr {other.m_end_data_ptr},
    m_num_bytes {other.m_num_bytes},
    m_type {other.m_type},
    m_current_data_value {other.m_current_data_value},
    m_is_current_pointee_cached {other.m_is_current_pointee_cached}
  {}

  /**
   * @copydoc IndividualFieldIterator::IndividualFieldIterator(const IndividualFieldIterator&)
   */
  IndividualFieldValueIterator(IndividualFieldValueIterator&& other) noexcept :
    m_body {std::move(other.m_body)},
    m_current_data_ptr {other.m_current_data_ptr},
    m_original_data_ptr {other.m_original_data_ptr},
    m_end_data_ptr {other.m_end_data_ptr},
    m_num_bytes {other.m_num_bytes},
    m_type {other.m_type},
    m_current_data_value {other.m_current_data_value},
    m_is_current_pointee_cached {other.m_is_current_pointee_cached}
  {}

  /**
   * @copydoc IndividualFieldIterator::IndividualFieldIterator(const IndividualFieldIterator&)
   */
  IndividualFieldValueIterator& operator=(const IndividualFieldValueIterator& other) {
    if (&this == other)
      return *this;
    m_body = std::move(other.m_body);
    m_current_data_ptr = other.m_current_data_ptr;
    m_original_data_ptr = other.m_original_data_ptr;
    m_end_data_ptr  = other.m_end_data_ptr;
    m_num_bytes = other.m_num_bytes;
    m_type = other.m_type;
    m_current_data_value = other.m_current_data_value;
    m_is_current_pointee_cached = other.m_is_current_pointee_cached;
    return *this;
  }


  /**
   * @copydoc IndividualFieldIterator::IndividualFieldIterator(IndividualFieldIterator&&)
   */
  IndividualFieldValueIterator& operator=(IndividualFieldValueIterator&& other) noexcept {
    if (&this == other)
      return *this;
    m_body = std::move(other.m_body);
    m_current_data_ptr = other.m_current_data_ptr;
    m_original_data_ptr = other.m_original_data_ptr;
    m_end_data_ptr  = other.m_end_data_ptr;
    m_num_bytes = other.m_num_bytes;
    m_type = other.m_type;
    m_current_data_value = other.m_current_data_value;
    m_is_current_pointee_cached = other.m_is_current_pointee_cached;
    return *this;
  }

  /**
   * @copydoc IndividualFieldIterator::operator+=(const int)
   */
  IndividualFieldValueIterator& operator+=(const int n) {
    m_current_data_ptr += n * m_num_bytes;
    m_is_current_pointee_cached = false;
    m_current_data_ptr = utils::cache_and_advance_to_end_if_necessary(m_current_data_ptr, m_end_data_ptr, *this);
    return *this;
  }

  /**
   * @copydoc IndividualFieldIterator::operator+=(int)
   */
  IndividualFieldValueIterator& operator-=(const int n) {
    m_current_data_ptr -= n * m_num_bytes;
    m_is_current_pointee_cached = false;
    return *this;
  }

  /**
   * @brief two iterators are equal if they are in exactly the same state (pointing at the same location in memory
   */
  bool operator==(const IndividualFieldValueIterator& other) { return (m_body == other.m_body && m_current_data_ptr == other.m_current_data_ptr); }

  /**
   * @brief the oposite check of IndividualFieldValueIterator::operator==()
   */
  bool operator!=(const IndividualFieldValueIterator& other) {
    return !(*this == other);
  }

  /**
   * @brief an operator is greater/less than another iterator if it is pointing to a previous element (sample) in the FormatField 
   * object. The order is determined by the Variant record.
   */
  bool operator<(const IndividualFieldValueIterator& other) {
    return m_body == other.m_body && m_current_data_ptr < other.m_current_data_ptr;
  }

  /**
   * @copydoc IndividualFieldValueIterator::operator<()
   */
  bool operator>(const IndividualFieldValueIterator& other) {
    return m_body == other.m_body && m_current_data_ptr > other.m_current_data_ptr;
  }

  /**
   * @copydoc IndividualFieldValueIterator::operator<()
   */
  bool operator<=(const IndividualFieldValueIterator& other) {
    return m_body == other.m_body && m_current_data_ptr <= other.m_current_data_ptr;
  }

  /**
   * @copydoc IndividualFieldValueIterator::operator<()
   */
  bool operator>=(const IndividualFieldValueIterator& other) {
    return m_body == other.m_body && m_current_data_ptr >= other.m_current_data_ptr;
  }

  /**
   * @brief direct access to the value of the current sample
   * @return the value in its native type
   */
  VALUE_TYPE operator*() const noexcept {
    if(m_is_current_pointee_cached)
      return m_current_data_value;
    else
      return convert_from_byte_array(m_current_data_ptr, 0);
  }

  /*
   * @brief - operator* is const function, hence cannot set m_current_data_value and m_is_current_pointee_cached
   */
  VALUE_TYPE read_and_cache_current_pointee() noexcept {
    m_current_data_value = *(*this);
    m_is_current_pointee_cached = true;
    return m_current_data_value;
  }
  /**
   * @brief advances to the next sample
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result of operator* with end(). This is the STL way.
   * @return modified iterator
   */
  IndividualFieldValueIterator& operator++() noexcept {
    m_current_data_ptr += m_num_bytes;
    m_is_current_pointee_cached = false;
    m_current_data_ptr = utils::cache_and_advance_to_end_if_necessary(m_current_data_ptr, m_end_data_ptr, *this);
    return *this;
  }

  /**
   * @brief postfix operator to advance iterator to the next position
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result of operator* with end(). This is the STL way.
   * @return copy of iterator at current position (before the increment)
   */
  IndividualFieldValueIterator operator++(int) noexcept {
    auto const tmp = IndividualFieldValueIterator(*this);
    operator++();
    return tmp;
  }
  /**
   * @brief advances to the previous sample
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result of operator* with end(). This is the STL way.
   * @return modified iterator
   */
  IndividualFieldValueIterator& operator--() {
    m_current_data_ptr -= m_num_bytes;
    m_is_current_pointee_cached = false;
    return *this;
  }

  /**
   * @brief postfix operator to retreat iterator to the previous position
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result of operator* with end(). This is the STL way.
   * @return copy of iterator at current position (before the decrement)
   */
  IndividualFieldValueIterator operator--(int) noexcept {
    auto const tmp = IndividualFieldValueIterator(*this);
    operator--();
    return tmp;
  }
  /**
   * @brief random access to the value of a given index for reading or writing
   * @param index must be between 0 and the number of indices for this record but no boundary check is done in this implementation
   * @note implementation guarantees this operation to be O(1)
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result of operator* with end(). This is the STL way.
   * @return the value in it's native type
   */
  VALUE_TYPE operator[](const uint32_t index) const {
    //if index == 0, use the cached value logic in operator*
    return (index == 0 ? *(*this) : convert_from_byte_array(m_current_data_ptr, index));
  }

 private:
  std::shared_ptr<bcf1_t> m_body;
  const uint8_t* m_current_data_ptr;
  const uint8_t* m_original_data_ptr;
  const uint8_t* m_end_data_ptr;
  uint32_t m_num_bytes;
  utils::VariantFieldType m_type;
  VALUE_TYPE m_current_data_value;
  bool m_is_current_pointee_cached;

  VALUE_TYPE convert_from_byte_array(const uint8_t* data_ptr, int index) const;
};

/**
 * @brief specialization of the conversion template for int32_t
 */
template<> inline
int32_t IndividualFieldValueIterator<int32_t>::convert_from_byte_array(const uint8_t* data_ptr, int index) const {
  return utils::convert_data_to_integer(data_ptr, index, m_num_bytes, m_type);
}

/**
 * @brief specialization of the conversion template for floats
 */
template<> inline
float IndividualFieldValueIterator<float>::convert_from_byte_array(const uint8_t* data_ptr, int index) const {
  return utils::convert_data_to_float(data_ptr, index, m_num_bytes, m_type);
}

/**
 * @brief specialization of the conversion template for strings
 */
template<> inline
std::string IndividualFieldValueIterator<std::string>::convert_from_byte_array(const uint8_t* data_ptr, int index) const {
  return utils::convert_data_to_string(data_ptr, index, m_num_bytes, m_type);
}

};

#endif // gamgee__individual_field_value_iterator__guard
