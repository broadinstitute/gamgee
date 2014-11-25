#ifndef gamgee__variant_builder_multi_sample_vector__guard
#define gamgee__variant_builder_multi_sample_vector__guard

#include <vector>
#include <string>
#include <cassert>

namespace gamgee {

/**
 * @brief Class that allows you to efficiently prepare multi-sample data for setting individual fields in VariantBuilder
 *
 * For those who want higher performance than is possible with the vector<vector> individual field setters
 * in VariantBuilder, this class is available as an alternative. Internally it uses a flattened, pre-padded
 * one-dimensional vector that can be passed directly to htslib for encoding, and so has much better
 * data locality and cache performance than a two-dimensional vector.
 *
 * To use, first determine the number of samples and the maximum number of values per
 * sample for the field. Then get a pre-initialized VariantBuilderMultiSampleVector
 * from your VariantBuilder instance. Eg.,
 *
 * auto multi_sample_vector = builder.get_integer_multi_sample_vector(num_samples, max_values_per_sample);
 *
 * (NOTE: DO NOT INSTANTIATE THIS CLASS DIRECTLY! ALWAYS GET AN INSTANCE FROM A VARIANTBUILDER OBJECT)
 *
 * This multi_sample_vector will have missing values for all samples, with appropriate padding to the
 * maximum field width.
 *
 * Then, fill in the values for each non-missing sample by invoking the set_sample_value() and/or
 * set_sample_values() functions on your multi-sample vector (NOTE: set_sample_value() is MUCH more efficient
 * than set_sample_values() since it doesn't require a vector construction/destruction for each call).
 * You don't have to worry about samples with no values, since all samples start out with missing values.
 *
 * Finally, pass your multi-sample vector into a VariantBuilder setter function (favoring the functions
 * that take field indices and use move semantics for high performance). Eg.,
 *
 * builder.set_integer_individual_field(field_index, std::move(multi_sample_vector));
 */
template<class ELEMENT_TYPE>
class VariantBuilderMultiSampleVector {
 public:

  /**
   * @brief Construct a VariantBuilderMultiSampleVector
   *
   * @param num_samples number of samples we will be storing in this vector
   * @param max_values_per_sample maximum number of values across all samples
   * @param missing_value the BCF missing value for this type
   * @param end_of_vector_value the BCF end of vector value for this type
   *
   * @note: DO NOT INSTANTIATE THIS CLASS DIRECTLY! ALWAYS GET AN INSTANCE FROM A VARIANTBUILDER OBJECT
   */
  VariantBuilderMultiSampleVector(const uint32_t num_samples, const uint32_t max_values_per_sample, const ELEMENT_TYPE missing_value, const ELEMENT_TYPE end_of_vector_value) :
    m_multi_sample_values(num_samples * max_values_per_sample, end_of_vector_value),  // Fill with end_of_vector_value
    m_num_samples{num_samples},
    m_max_values_per_sample{max_values_per_sample}
  {
    // Place a single missing value at the start of each sample's values (the rest of the vector has already
    // been padded with vector end values)
    for ( auto sample_start = 0u; sample_start < m_multi_sample_values.size(); sample_start += m_max_values_per_sample ) {
      m_multi_sample_values[sample_start] = missing_value;
    }
  }

  // Both copyable and moveable, with default destruction
  VariantBuilderMultiSampleVector(const VariantBuilderMultiSampleVector& other) = default;
  VariantBuilderMultiSampleVector& operator=(const VariantBuilderMultiSampleVector& other) = default;
  VariantBuilderMultiSampleVector(VariantBuilderMultiSampleVector&& other) = default;
  VariantBuilderMultiSampleVector& operator=(VariantBuilderMultiSampleVector&& other) = default;
  ~VariantBuilderMultiSampleVector() = default;

  /**
   * @brief Set a single value for one sample
   *
   * @param sample_index index of the sample (from a VariantHeader lookup)
   * @param value_index index of the value we are setting for this sample (values for EACH sample start at index 0)
   * @param value value to set
   *
   * @note MUCH more efficient than set_sample_values() below, since it doesn't require a vector
   *       construction/destruction for each call.
   *
   * @warning Bounds checking is performed only in debug builds (for the sake of performance), so be sure that sample_index
   *          is < num_samples and value_index < max_values_per_sample
   */
  inline void set_sample_value(const uint32_t sample_index, const uint32_t value_index, const ELEMENT_TYPE value) {
    assert(sample_index < m_num_samples);
    assert(value_index < m_max_values_per_sample);

    m_multi_sample_values[sample_index * m_max_values_per_sample + value_index] = value;
  }

  /**
   * @brief Set all values for one sample at once
   *
   * @param sample_index index of the sample (from a VariantHeader lookup)
   * @param values vector of values for this sample
   *
   * @note LESS efficient than setting one value at a time using set_sample_value(), since this function
   *       involves creating/destroying a vector for each sample.
   *
   * @warning Bounds checking is performed only in debug builds (for the sake of performance), so be sure that sample_index
   *          is < num_samples and values.size() <= max_values_per_sample
   */
  inline void set_sample_values(const uint32_t sample_index, const std::vector<ELEMENT_TYPE>& values) {
    assert(sample_index < m_num_samples);
    assert(values.size() <= m_max_values_per_sample);

    const auto sample_start = sample_index * m_max_values_per_sample;
    for ( auto value_index = 0u; value_index < values.size(); ++value_index ) {
      m_multi_sample_values[sample_start + value_index] = values[value_index];
    }
  }

  /**
   * @brief Get a reference to the internal one-dimensional vector used for value storage
   */
  const std::vector<ELEMENT_TYPE>& get_vector() const { return m_multi_sample_values; }

  /**
   * @brief Return the number of samples whose values we are storing
   */
  uint32_t num_samples() const { return m_num_samples; }

  /**
   * @brief Return the maximum number of values across all samples
   */
  uint32_t max_values_per_sample() const { return m_max_values_per_sample; }

 private:
  std::vector<ELEMENT_TYPE> m_multi_sample_values;
  uint32_t m_num_samples;
  uint32_t m_max_values_per_sample;

  friend class VariantBuilder; // VariantBuilder needs access to internals in order to build efficiently
};

}

#endif  /* gamgee__variant_builder_multi_sample_vector__guard */