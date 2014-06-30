#ifndef gamgee__short_value_optimized_storage__guard
#define gamgee__short_value_optimized_storage__guard

#include "utils.h"

#include <vector>
#include <utility>
#include <cstring>
#include <stdexcept>

namespace gamgee {
namespace utils {

/**
 * A class for storing arrays of values of varying lengths, optimizing for short arrays (where "short" is custom-defined).
 * "Short" arrays are stored in preallocated contiguous storage, while arrays that don't fit in contiguous storage
 * are dynamically allocated. The idea is to preserve data locality and avoid new memory allocations/deallocations as
 * much as possible. Performance should be good provided that most arrays fall into the "short" category.
 *
 * This is intended as a more-efficient replacement for vector<vector> in cases where typical inner vector
 * lengths are known in advance.
 */
template<class ELEMENT_TYPE>
class ShortValueOptimizedStorage {
 public:
  ShortValueOptimizedStorage(const uint32_t capacity, const uint32_t short_value_upper_bound) :
    m_num_values { 0 },
    m_max_value_length { 0 },
    m_short_value_upper_bound { short_value_upper_bound },
    m_values(capacity),
    m_contiguous_storage(capacity * short_value_upper_bound)
  {
    if ( short_value_upper_bound == 0 ) {
      throw std::invalid_argument{"short value upper bound must be > 0"};
    }

    for ( auto i = 0u; i < m_values.size(); ++i ) {
      m_values[i].first = nullptr;
      m_values[i].second = 0u;
    }
  }

  ~ShortValueOptimizedStorage() {
    // Custom destructor because we have conditional deletion of values in the container
    clear();
  }

  // Moveable but not copyable
  ShortValueOptimizedStorage(ShortValueOptimizedStorage&& other) = default;
  ShortValueOptimizedStorage& operator=(ShortValueOptimizedStorage&& other) = default;
  ShortValueOptimizedStorage(const ShortValueOptimizedStorage& other) = delete;
  ShortValueOptimizedStorage& operator=(const ShortValueOptimizedStorage& other) = delete;

  /**
   * @brief Returns the number of values in the container
   */
  uint32_t num_values() const { return m_num_values; }

  /**
   * @brief Returns the maximum number of values the container can hold
   */
  uint32_t capacity() const { return m_values.size(); }

  /**
   * @brief Returns the length of the longest value in the container
   */
  uint32_t max_value_length() const { return m_max_value_length; }

  /**
   * @brief Returns the length of the value at the specified position
   */
  uint32_t value_length(const uint32_t index) const {
    check_max_boundary(index, capacity());
    return m_values[index].second;
  }

  /**
   * @brief Determines whether there is a value at the specified position
   */
  bool is_set(const uint32_t index) const {
    return index < capacity() && m_values[index].first != nullptr;
  }

  /**
   * @brief Gets the value stored at a specified index (first of pair) plus its length in number of items (second of pair)
   */
  std::pair<ELEMENT_TYPE*,uint32_t> get(const uint32_t index) const {
    check_max_boundary(index, capacity());
    return m_values[index];
  }

  /**
   * @brief Set the value at the specified index by vector
   */
  void set(const uint32_t index, const std::vector<ELEMENT_TYPE>& values) {
    if ( values.empty() ) return;
    set(index, &(values[0]), values.size());
  }

  /**
   * @brief Set the value at the specified index by pointer
   */
  void set(const uint32_t index, const ELEMENT_TYPE* values, const uint32_t num_values) {
    // Do nothing for empty values
    if ( num_values == 0 || values == nullptr ) return;

    // Make sure index is in bounds
    check_max_boundary(index, capacity());

    // If this is the first time we've used this slot, we have one more value in the container
    if ( m_values[index].first == nullptr ) { ++m_num_values; }

    // Determine how the value will be stored based on its length, and get storage ready
    prepare_memory(index, num_values);

    // Copy values into the prepared storage
    memcpy(m_values[index].first, values, num_values * sizeof(ELEMENT_TYPE));

    // Keep track of the length of the longest set of values in the container
    if ( num_values > m_max_value_length) m_max_value_length = num_values;
  }

  /**
   * @brief Clear the value at a specific index
   *
   * @note may require us to recalculate the length of the longest value in storage
   */
  void clear(const uint32_t index) {
    if ( ! is_set(index) ) return;

    const auto previous_length = m_values[index].second;
    release_memory_if_necessary(index);
    --m_num_values;

    // Special case: clearing the longest item requires us to recalculate m_max_value_length
    if ( previous_length == m_max_value_length ) {
      recalculate_max_value_length();
    }
  }

  /**
   * @brief Reset storage to a pristine state with no values
   */
  void clear() {
    for ( auto i = 0u; i < m_values.size(); ++i ) {
      release_memory_if_necessary(i);
    }
    m_num_values = 0;
    m_max_value_length = 0;
  }

 private:
  uint32_t m_num_values;                                     ///< number of values in the storage
  uint32_t m_max_value_length;                               ///< length of the longest value in the storage
  uint32_t m_short_value_upper_bound;                        ///< length of the longest "short" value; values longer than this will be stored dynamically
  std::vector<std::pair<ELEMENT_TYPE*,uint32_t>> m_values;   ///< pointer to data + data length, for each field; cannot use smart pointers here, since we may either be pointing into contiguous storage or dynamically-allocated storage
  std::vector<ELEMENT_TYPE> m_contiguous_storage;            ///< pre-allocated storage used for short values; m_values may point into here

  /**
   * @brief Prepares storage for a new value by choosing between contiguous storage and dynamic allocation
   *        based on the value length
   */
  void prepare_memory(const uint32_t index, const uint32_t num_values) {
    // Clean up previous dynamic allocation at this index if there is one
    release_memory_if_necessary(index);

    // Use contiguous pre-allocated storage for short values
    if ( num_values <= m_short_value_upper_bound ) {
      m_values[index].first = &(m_contiguous_storage[index * m_short_value_upper_bound]);
    }
    else {
      // Use dynamic allocation for values that don't fit (note: can't use smart pointers to wrap
      // this allocation since we will often be pointing into memory that must not be freed -- otherwise
      // we risk undefined behavior)
      m_values[index].first = new ELEMENT_TYPE[num_values];
    }

    // Record the length of the value as well
    m_values[index].second = num_values;
  }

  /**
   * @brief Frees allocated memory at an index, but only if it was dynamically allocated instead of
   *        pointing into contiguous storage
   */
  void release_memory_if_necessary(const uint32_t index) {
    // Only free pointers that don't point into the contiguous storage
    if ( m_values[index].first != nullptr &&
         m_values[index].first != &(m_contiguous_storage[index * m_short_value_upper_bound]) ) {

      // No complaining about this delete! If you're the author of a low-level container like vector,
      // or this class, delete is perfectly ok to use internally. Can't use smart pointers because we
      // have pointers into contiguous storage mixed in with pointers to individually-allocated memory
      // blocks.
      delete[] m_values[index].first;
    }

    m_values[index].first = nullptr;
    m_values[index].second = 0;
  }

  /**
   * @brief Occasionally (ie., after removals) we need to examine all values to re-calculate
   *        the length of the longest value
   */
  void recalculate_max_value_length() {
    m_max_value_length = 0u;

    for ( auto& value : m_values ) {
      if ( value.first != nullptr && value.second > m_max_value_length ) {
        m_max_value_length = value.second;
      }
    }
  }
};

}
}

#endif  /* gamgee__short_value_optimized_storage__guard */
