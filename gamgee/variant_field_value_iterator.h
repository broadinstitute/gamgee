#ifndef gamgee__format_field_generic_value_iterator__
#define gamgee__format_field_generic_value_iterator__

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
  /**
   * @copydoc FormatFieldIterator::FormatFieldIterator(const FormatFieldIterator&)
   */
  VariantFieldValueIterator& operator=(const VariantFieldValueIterator& other) = delete;

  /**
   * @copydoc FormatFieldIterator::FormatFieldIterator(FormatFieldIterator&&)
   */
  VariantFieldValueIterator& operator=(VariantFieldValueIterator&& other) noexcept {
    if (&this == other)
      return *this;
    m_body = std::move(other.m_body);
    m_data_ptr = other.m_data_ptr;
  }

  /**
   * @copydoc FormatFieldIterator::operator+=(const int)
   */
  VariantFieldValueIterator& operator+=(const int n) {
    m_data_ptr += n;
    return *this;
  }

  /**
   * @copydoc FormatFieldIterator::operator+=(int)
   */
  VariantFieldValueIterator& operator-=(const int n) {
    m_data_ptr -= n;
    return *this;
  }

  /**
   * @brief two iterators are equal if they are in exactly the same state (pointing at the same location in memory
   */
  bool operator==(const VariantFieldValueIterator& other) {
    return m_body == other.m_body && m_data_ptr == other.m_data_ptr;
  }

  /**
   * @brief the oposite check of VariantFieldValueIterator::operator==()
   */
  bool operator!=(const VariantFieldValueIterator& other) {
    return m_body != other.m_body || m_data_ptr != other.m_data_ptr;
  }

  /**
   * @brief an operator is greater/less than another iterator if it is pointing to a previous element (sample) in the FormatField 
   * object. The order is determined by the Variant record.
   */
  bool operator<(const VariantFieldValueIterator& other) {
    return m_body == other.m_body && m_data_ptr < other.m_data_ptr;
  }

  /**
   * @copydoc VariantFieldValueIterator::operator<()
   */
  bool operator>(const VariantFieldValueIterator& other) {
    return m_body == other.m_body && m_data_ptr > other.m_data_ptr;
  }

  /**
   * @copydoc VariantFieldValueIterator::operator<()
   */
  bool operator<=(const VariantFieldValueIterator& other) {
    return m_body == other.m_body && m_data_ptr <= other.m_data_ptr;
  }

  /**
   * @copydoc VariantFieldValueIterator::operator<()
   */
  bool operator>=(const VariantFieldValueIterator& other) {
    return m_body == other.m_body && m_data_ptr >= other.m_data_ptr;
  }

  /**
   * @brief direct access to the value of the current sample
   * @return the value in its native type
   */
  VALUE_TYPE& operator*() const noexcept {
    return VALUE_TYPE{m_body, m_data_ptr};
  }

  /**
   * @brief advances to the next sample
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result of operator* with end(). This is the STL way.
   * @return the next value in it's native type
   */
  VALUE_TYPE& operator++() noexcept {
    ++m_data_ptr;
    return *m_data_ptr;
  }

  /**
   * @brief advances to the previous sample
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result of operator* with end(). This is the STL way.
   * @return the previous value in it's native type
   */
  VALUE_TYPE& operator--() {
    --m_data_ptr;
    return *m_data_ptr;
  }

  /**
   * @brief random access to the value of a given index for reading or writing
   * @param index must be between 0 and the number of indices for this record but no boundary check is done in this implementation
   * @note implementation guarantees this operation to be O(1)
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result of operator* with end(). This is the STL way.
   * @return the value in it's native type
   */
  VALUE_TYPE& operator[](const uint32_t index) const {
    return *m_data_ptr[index];
  }

 private:
  const std::shared_ptr<bcf1_t> m_body;
  VALUE_TYPE* m_data_ptr;
};

}

#endif
