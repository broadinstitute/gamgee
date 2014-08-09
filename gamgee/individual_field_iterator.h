#ifndef gamgee__individual_field_iterator__guard
#define gamgee__individual_field_iterator__guard 

#include "utils/utils.h"

#include "htslib/vcf.h"

#include<iterator>

namespace gamgee {


/** 
 * @brief iterator for VariantField objects. 
 * 
 * This iterator will walk through all the samples in a Variant record for a
 * given VariantField object. For example if you want to iterate over all the GQ
 * values of a Variant record you would do so through this iterator.
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
template<class TYPE>
class IndividualFieldIterator : public std::iterator<std::random_access_iterator_tag, TYPE> {
 public:

  /**
   * @brief simple constructor used by VariantField to create an iterator
   * @param body the Variant record used to access the data (shared ownership guarantees availability until iterator goes out of scope)
   * @param format_ptr pointer to the raw byte array in body where all the values for this format field is encoded
   * @param end_iterator whether or not this is being called by the VariantField::end() function.
   * @note this constructor serves only the VariantField::begin() and VariantField::end() functions. 
   */
  IndividualFieldIterator(const std::shared_ptr<bcf1_t>& body, const bcf_fmt_t* const format_ptr, bool end_iterator = false) :
    m_body {body}, 
    m_format_ptr {format_ptr},
    m_data_ptr {end_iterator ? format_ptr->p + m_format_ptr->size * m_body->n_sample : format_ptr->p} 
  {}

  /**
   * @brief copy constructor is only meant for internal STL functions use. It makes a shallow copy
   * of the underlying object which is sufficient for an iterator, but not exactly what a user would
   * expect (user should expect copy constructors to make deep copies). But a deep copied iterator
   * makes no sense.
   * @warning does not deep copy the underlying data! (but the copied iterator will be able to
   * navigate independentely). 
   */
  IndividualFieldIterator(const IndividualFieldIterator& other) :
    m_body {other.m_body},
    m_format_ptr {other.m_format_ptr},
    m_data_ptr {other.m_data_ptr}
  {}
      
  /**
   * @brief safely moves the data from one VariantField to a new one without making any copies
   * @param other another VariantField object
   */
  IndividualFieldIterator(IndividualFieldIterator&& other) noexcept : 
    m_body {std::move(other.m_body)},
    m_format_ptr {other.m_format_ptr},
    m_data_ptr {other.m_data_ptr}
  {}

  /**
   * @copydoc IndividualFieldIterator::IndividualFieldIterator(const IndividualFieldIterator&)
   */
  IndividualFieldIterator& operator=(const IndividualFieldIterator& other) {
    if (&this == other)
      return *this;
    m_body = other.m_body;
    m_format_ptr = other.m_format_ptr;
    m_data_ptr = other.m_data_ptr;
    return *this;
  }

  /**
   * @copydoc IndividualFieldIterator::IndividualFieldIterator(IndividualFieldIterator&&)
   */
  IndividualFieldIterator& operator=(IndividualFieldIterator&& other) noexcept {
    if (&this == other)
      return *this;
    m_body = std::move(other.m_body);
    m_format_ptr = other.m_format_ptr;
    m_data_ptr = other.m_data_ptr;
    return *this;
  }

  /**
   * @brief simple compound assignment operation for random advances (back/forward) to the iterator
   * @param n how much to advance (negative numbers to go the other direction)
   * @warning there is no boundary check in this operator
   */
  IndividualFieldIterator& operator+=(const int n) noexcept {
    m_data_ptr += n * m_format_ptr->size;
    return *this;
  }

  /**
   * @copydoc IndividualFieldIterator::operator+=(int)
   */
  IndividualFieldIterator& operator-=(const int n) noexcept {
    m_data_ptr -= n * m_format_ptr->size;
    return *this;
  }

  /**
   * @brief two iterators are equal if they are in exactly the same state (pointing at the same location in memory
   */
  bool operator==(const IndividualFieldIterator& other) {
    return m_body == other.m_body && m_data_ptr == other.m_data_ptr;
  }

  /**
   * @brief the oposite check of IndividualFieldIterator::operator==()
   */
  bool operator!=(const IndividualFieldIterator& other) {
    return m_body != other.m_body || m_data_ptr != other.m_data_ptr;
  }

  /**
   * @brief an operator is greater/less than another iterator if it is pointing to a previous element (sample) in the VariantField 
   * object. The order is determined by the Variant record.
   */
  bool operator<(const IndividualFieldIterator& other) {
    return m_body == other.m_body && m_data_ptr < other.m_data_ptr;
  }

  /**
   * @copydoc IndividualFieldIterator::operator<()
   */
  bool operator>(const IndividualFieldIterator& other) {
    return m_body == other.m_body && m_data_ptr > other.m_data_ptr;
  }

  /**
   * @copydoc IndividualFieldIterator::operator<()
   */
  bool operator<=(const IndividualFieldIterator& other) {
    return m_body == other.m_body && m_data_ptr <= other.m_data_ptr;
  }

  /**
   * @copydoc IndividualFieldIterator::operator<()
   */
  bool operator>=(const IndividualFieldIterator& other) {
    return m_body == other.m_body && m_data_ptr >= other.m_data_ptr;
  }

  /**
   * @brief direct access to the value of the current sample
   * @return the value if it is a basic type (e.g. GQ, GL), or a specific object if it is a complex type (e.g. PL, AD,...)
   */
  TYPE operator*() const noexcept {
    return TYPE{m_body, m_format_ptr, m_data_ptr};
  }

  /**
   * @brief Prefix increment. Advances to the next sample
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result with end(). This is the STL way.
   * @return a reference to the start of the values of the next sample
   */
  IndividualFieldIterator& operator++() noexcept {
    operator+=(1);
    return *this;
  }

  /**
   * @brief Postfix increment. Advances to the next sample
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the end by comparing the result with end(). This is the STL way.
   * @return the value of the start of the same sample.
   */
  IndividualFieldIterator operator++(int) noexcept {
    auto const tmp = IndividualFieldIterator(*this);
    operator++();
    return tmp;
  }

  /**
   * @brief Prefix increment.  Reverses to the previous sample.
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the beginning by comparing the result with begin(). This is the STL way.
   * @return a reference to the start of the values of the previous sample
   */
  IndividualFieldIterator& operator--() noexcept {
    operator-=(1);
    return *this;
  }

  /**
   * @brief Postfix increment.  Reverses to the previous sample.
   * @note mainly designed for iterators
   * @warning does not check for bounds exception, you should verify whether or not you've reached the beginning by comparing the result with begin(). This is the STL way.
   * @return the value of the start of the same sample.
   */
  IndividualFieldIterator& operator--(int) noexcept {
    auto const tmp = IndividualFieldIterator(*this);
    operator--();
    return tmp;
  }

  /**
   * @brief random access to the value of a given sample for reading or writing
   * @param sample must be between 0 and the number of samples for this record 
   * @note implementation guarantees this operation to be O(1)
   * @exception std::out_of_range if index is out of range
   * @return the value if it is a basic type (e.g. GQ, GL), or a specific object if it is a complex type (e.g. PL, AD,...)
   */
  TYPE operator[](const uint32_t sample) const {
    utils::check_max_boundary(sample, m_body->n_sample);
    return TYPE{m_body, m_format_ptr, m_format_ptr->p + (sample * m_format_ptr->size)};
  }

  /**
   * @brief Difference between two iterators as an integer.
   * @note Useful as a substitute for size() when only begin() and end() are available.
   * @returns the number of iterator steps between [first, last) where last is the current IndividualFieldIterator.
   * @param first is the iterator the position of which is to be subtracted from the position of the current iterator.
   */
  int32_t  operator-(const IndividualFieldIterator<TYPE>& first) const {
    return static_cast<int32_t>(m_data_ptr - first.m_data_ptr)/m_format_ptr->size;
  }

 private:
  std::shared_ptr<bcf1_t> m_body;      ///< shared ownership of the Variant record memory so it stays alive while this object is in scope
  const bcf_fmt_t* const m_format_ptr; ///< pointer to the format_field in the body so we can access the tag's information
  uint8_t* m_data_ptr;                 ///< pointer to m_body structure where the data for this particular type is located.

};

}

#endif // gamgee__individual_field_iterator__guard
