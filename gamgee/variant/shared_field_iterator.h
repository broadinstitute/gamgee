#ifndef gamgee__shared_field_iterator__guard
#define gamgee__shared_field_iterator__guard

#include "../utils/variant_field_type.h"

#include "htslib/vcf.h"

#include <iterator>
#include <memory>

namespace gamgee {

/** 
 * @brief iterator for SharedField objects. 
 * 
 * This iterator will walk through all the values for a given
 * SharedField object. For example if you want to iterate over all
 * the AN values of a Variant record you would do so through this iterator.
 *
 * @note implements a random access iterator which gives you full performance
 * on STL algorithms that use iterators (mostly every one)
 * 
 * @note this iterator never makes any copies of the underlying memory. It is
 * implemented with pointer acrobatics behind the scenes for maximum
 * performance while maintaining a friendly interface.
 */
template<class VALUE_TYPE>
class SharedFieldIterator : public std::iterator<std::random_access_iterator_tag, VALUE_TYPE> {
 public:

  /**
   * @brief default constructor of an empty iterator
   * @note useful for creating begin/end iterators over empty objects that match
   */
  SharedFieldIterator() : 
    m_body {nullptr},
    m_current_data_ptr {nullptr},
    m_original_data_ptr{nullptr},
    m_end_data_ptr{nullptr},
    m_bytes_per_value{0},
    m_type{utils::VariantFieldType::NIL},
    m_is_current_pointee_cached{false}
  {}
 
  /**
   * @brief Constructor with bcf1_t structure and start and end pointers of the vector/array
   * @param body a pointer to the bcf1_t data structure to be held as a shared pointer
   * @param data_ptr the byte array containing the value(s) for this shared field
   * @param end_ptr the end pointer value of the data array
   * @param bytes_per_value number of bytes in each value
   * @param type the encoding of the value 
   * @note this constructor is probably only used by SharedField::begin() and SharedField::end()
   */  
  explicit SharedFieldIterator(const std::shared_ptr<bcf1_t>& body, uint8_t* data_ptr, uint8_t* end_ptr, const uint32_t bytes_per_value, const utils::VariantFieldType& type) :
    m_body {body},
    m_current_data_ptr {data_ptr},
    m_original_data_ptr {data_ptr},
    m_end_data_ptr {end_ptr},
    m_bytes_per_value {bytes_per_value},
    m_type {type},
    m_is_current_pointee_cached{false}
  {}

  /**
   * @brief Constructor with bcf1_t structure and start pointer of the vector/array
   * @param body a pointer to the bcf1_t data structure to be held as a shared pointer
   * @param data_ptr the byte array containing the value(s) for this shared field
   * @param bytes_per_value number of bytes in each value
   * @param type the encoding of the value 
   * @note this constructor is probably only used by SharedField::begin() and SharedField::end()
   */
  explicit SharedFieldIterator(const std::shared_ptr<bcf1_t>& body, uint8_t* data_ptr, const uint32_t bytes_per_value, const utils::VariantFieldType& type): SharedFieldIterator(body, data_ptr, nullptr, bytes_per_value, type)
  {}
   
  SharedFieldIterator(const SharedFieldIterator& other) = default; ///< standard copy constructor creates a new iterator pointing to the same underlying data
  SharedFieldIterator(SharedFieldIterator&& other) = default; ///< standard move constructor
  SharedFieldIterator& operator=(const SharedFieldIterator& other) = default;  ///< standard copy assignment operator creates a new iterator pointing to the same underlying data
  SharedFieldIterator& operator=(SharedFieldIterator&& other) = default; ///< standard move assignment operator
  ///< two iterators are equal if they are in exactly the same state (pointing at the same location in memory
  bool operator==(const SharedFieldIterator& other) { return (m_body == other.m_body && m_current_data_ptr == other.m_current_data_ptr); } 
  bool operator!=(const SharedFieldIterator& other) { return !(*this == other); }                                                       ///< the negation of SharedFieldIterator::operator==()
  bool operator<(const SharedFieldIterator& other) { return m_body == other.m_body && m_current_data_ptr < other.m_current_data_ptr;}   ///< an operator is greater/less than another iterator if it is pointing to a previous element in the SharedField object. The order is determined by the Variant record.
  bool operator>(const SharedFieldIterator& other) { return m_body == other.m_body && m_current_data_ptr > other.m_current_data_ptr;}    ///< not smaller than other neither equal to other
  bool operator<=(const SharedFieldIterator& other) { return m_body == other.m_body && m_current_data_ptr <= other.m_current_data_ptr;} ///< not greater than other
  bool operator>=(const SharedFieldIterator& other) { return m_body == other.m_body && m_current_data_ptr >= other.m_current_data_ptr;} ///< not smaller than other

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

  /** advance n values */
  SharedFieldIterator& operator+=(const int n) {
    m_current_data_ptr += n * m_bytes_per_value;
    m_is_current_pointee_cached = false;
    m_current_data_ptr = utils::cache_and_advance_to_end_if_necessary(m_current_data_ptr, m_end_data_ptr, *this);
    return *this;
  }

  /** moves back n values */
  SharedFieldIterator& operator-=(const int n) {
    m_current_data_ptr -= n * m_bytes_per_value;
    m_is_current_pointee_cached = false;
    return *this;
  }

  /**
   * @brief advances to the next value
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result of operator* with end(). This is the STL way.
   * @return updated iterator
   */
  SharedFieldIterator& operator++() noexcept {
    m_current_data_ptr += m_bytes_per_value;
    m_is_current_pointee_cached = false;
    m_current_data_ptr = utils::cache_and_advance_to_end_if_necessary(m_current_data_ptr, m_end_data_ptr, *this);
    return *this;
  }

  /**
   * @brief Postfix increment. Advances to the next sample
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result with end(). This is the STL way.
   * @return copy of current iterator before the increment
   */
  SharedFieldIterator operator++(int) noexcept {
    const auto tmp = SharedFieldIterator(*this);
    operator++();
    return tmp;
  }


  /**
   * @brief advances to the previous value
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result of operator* with end(). This is the STL way.
   * @return updated iterator
   */
  SharedFieldIterator& operator--() {
    m_current_data_ptr -= m_bytes_per_value;
    m_is_current_pointee_cached = false;
    return *this;
  }

  /**
   * @brief Postfix decrement. Move back to the previous sample
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result with end(). This is the STL way.
   * @return copy of current iterator before the decrement
   */
  SharedFieldIterator operator--(int) noexcept {
    const auto tmp = SharedFieldIterator(*this);
    operator--();
    return tmp;
  }
  /**
   * @brief difference between two iterators as an integer.
   * @param first is the iterator the position of which is to be subtracted from the position of the current iterator.
   * @returns the number of values steps between first and last.
   */
  int32_t operator-(const SharedFieldIterator& first) const {
    return static_cast<int32_t>(m_current_data_ptr - first.m_current_data_ptr)/m_bytes_per_value;
  }

  /**
   * @brief random access to the value of a given index for reading or writing
   * @param index must be between 0 and the number of indices for this record but no boundary check is done in this implementation
   * @note implementation guarantees this operation to be O(1)
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result of operator* with end(). This is the STL way (for iterators).
   * @return the value in it's native type
   */
  VALUE_TYPE operator[](const uint32_t index) const {
    return (index == 0 ? *(*this) : convert_from_byte_array(m_current_data_ptr, index));
  }

 private:
  std::shared_ptr<bcf1_t> m_body;
  const uint8_t* m_current_data_ptr;
  const uint8_t* m_original_data_ptr;
  const uint8_t* m_end_data_ptr;
  uint32_t m_bytes_per_value;
  utils::VariantFieldType m_type;
  VALUE_TYPE m_current_data_value;
  bool m_is_current_pointee_cached;

  VALUE_TYPE convert_from_byte_array(const uint8_t* data_ptr, int index) const;
};

/**
 * @brief specialization of the conversion template for int32_t
 */
template<> inline
int32_t SharedFieldIterator<int32_t>::convert_from_byte_array(const uint8_t* data_ptr, int index) const {
  return utils::convert_data_to_integer(data_ptr, index, m_bytes_per_value, m_type);
}

/**
 * @brief specialization of the conversion template for floats
 */
template<> inline
float SharedFieldIterator<float>::convert_from_byte_array(const uint8_t* data_ptr, int index) const {
  return utils::convert_data_to_float(data_ptr, index, m_bytes_per_value, m_type);
}

/**
 * @brief specialization of the conversion template for strings
 */
template<> inline
std::string SharedFieldIterator<std::string>::convert_from_byte_array(const uint8_t* data_ptr, int index) const {
  return utils::convert_data_to_string(data_ptr, index, m_bytes_per_value, m_type);
}

};

#endif // gamgee__shared_field_iterator__guard
