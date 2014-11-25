#ifndef gamgee__variant_builder_individual_field__guard
#define gamgee__variant_builder_individual_field__guard

#include "htslib/kstring.h"
#include "htslib/vcf.h"

#include "../missing.h"
#include "../utils/hts_memory.h"
#include "../utils/short_value_optimized_storage.h"

#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace gamgee {

/**
 * @brief Helper class for VariantBuilder to manage the storage and encoding of a single multi-sample individual field
 *
 * Prior to build time, field data is stored in un-encoded form in one of two ways:
 *
 * -in vectors moved or copied from the user (these represent bulk changes to the entire field)
 *
 * -in ShortValueOptimizedStorage when there are per-sample changes -- this allows us to efficiently
 *  handle the setting of many small per-sample values without performing any extra dynamic memory
 *  allocation most of the time
 *
 * There cannot be both bulk and per-sample changes -- this is treated as an error, since it would be far
 * too expensive to reconcile the two.
 *
 * At build time, the raw, un-encoded field data is passed to htslib for encoding into the final byte
 * array for insertion into a Variant object.
 */
template<class ENCODED_TYPE, class BULK_CHANGE_TYPE>
class VariantBuilderIndividualField {
 public:
  explicit VariantBuilderIndividualField(const uint32_t num_samples, const uint32_t field_index, const int32_t field_type,
                                         const ENCODED_TYPE missing_value, const ENCODED_TYPE end_of_vector_value,
                                         const uint32_t short_value_upper_bound) :
    m_field_index {field_index},
    m_field_type {field_type},
    m_num_samples {num_samples},
    m_max_sample_value_length { 0 },
    m_missing_value {missing_value},
    m_end_of_vector_value {end_of_vector_value},
    m_flattened_bulk_changes{},
    m_nested_bulk_changes{},
    m_per_sample_changes {num_samples, short_value_upper_bound},
    m_removed { false }
  {}

  VariantBuilderIndividualField(VariantBuilderIndividualField&& other) = default;
  VariantBuilderIndividualField& operator=(VariantBuilderIndividualField&& other) = default;
  VariantBuilderIndividualField(const VariantBuilderIndividualField& other) = delete;
  VariantBuilderIndividualField& operator=(const VariantBuilderIndividualField& other) = delete;
  ~VariantBuilderIndividualField() = default;

  /**
   * Bulk setters, allowing field to be set for all samples using one- or
   * two-dimensional vectors, and by move or by copy
   */

  void set_entire_field(std::vector<BULK_CHANGE_TYPE>&& bulk_changes) {
    m_flattened_bulk_changes = std::move(bulk_changes);
    if ( ! m_nested_bulk_changes.empty() ) m_nested_bulk_changes.clear();
    m_max_sample_value_length = max_sample_value_length(m_flattened_bulk_changes);
    m_removed = false;
  }

  void set_entire_field(const std::vector<BULK_CHANGE_TYPE>& bulk_changes) {
    m_flattened_bulk_changes = bulk_changes;
    if ( ! m_nested_bulk_changes.empty() ) m_nested_bulk_changes.clear();
    m_max_sample_value_length = max_sample_value_length(m_flattened_bulk_changes);
    m_removed = false;
  }

  void set_entire_field(std::vector<std::vector<BULK_CHANGE_TYPE>>&& bulk_changes) {
    m_nested_bulk_changes = std::move(bulk_changes);
    if ( ! m_flattened_bulk_changes.empty() ) m_flattened_bulk_changes.clear();
    m_max_sample_value_length = max_sample_value_length(m_nested_bulk_changes);
    m_removed = false;
  }
  void set_entire_field(const std::vector<std::vector<BULK_CHANGE_TYPE>>& bulk_changes) {
    m_nested_bulk_changes = bulk_changes;
    if ( ! m_flattened_bulk_changes.empty() ) m_flattened_bulk_changes.clear();
    m_max_sample_value_length = max_sample_value_length(m_nested_bulk_changes);
    m_removed = false;
  }

  /**
   * @brief Stores a value for just a single sample efficiently using the ShortValueOptimizedStorage layer
   */
  void set_sample_field_value(const uint32_t sample_index, const ENCODED_TYPE* values, const uint32_t num_values) {
    m_per_sample_changes.set(sample_index, values, num_values);
    m_max_sample_value_length = m_per_sample_changes.max_value_length();
    m_removed = false;
  }

  uint32_t field_index() const { return m_field_index; }
  int32_t field_type() const { return m_field_type; }

  void remove() { clear(); m_removed = true; }
  bool removed() const { return m_removed; }
  bool missing() const { return m_max_sample_value_length == 0u; }
  bool present() const { return ! missing() && ! removed(); }

  bool has_bulk_changes() const { return ! m_flattened_bulk_changes.empty() || ! m_nested_bulk_changes.empty(); }
  bool has_per_sample_changes() const { return m_per_sample_changes.num_values() > 0; }

  /**
   * @brief Reset this field to a pristine state with no data
   */
  void clear() {
    if ( has_bulk_changes() ) {
      m_flattened_bulk_changes.clear();
      m_nested_bulk_changes.clear();
    }

    if ( has_per_sample_changes() ) {
      m_per_sample_changes.clear();
    }
    m_max_sample_value_length = 0;

    m_removed = false;
  }

  /**
   * @brief Provide an estimate (typically an overestimate) of the number of bytes this field will
   *        require when encoded
   */
  uint32_t estimated_encoded_size() const {
    // max possible size required by encoded field index + typing byte
    constexpr uint32_t max_metadata_overhead = 11;

    // note that for ints we (typically) overestimate by assuming the worst-case of int32
    return present() ? (m_max_sample_value_length * sizeof(ENCODED_TYPE) * m_num_samples) + max_metadata_overhead : 0;
  }

  /**
   * @brief Encode this field's data into the provided buffer. If field has no data or was removed, do nothing.
   */
  void encode_into(kstring_t* destination) const {
    if ( ! present() ) {
      return;
    }

    // Branch based on how the un-encoded field values are being stored
    if ( has_bulk_changes() && has_per_sample_changes() ) {
      throw std::logic_error("Cannot set an individual field both in bulk and by sample");
    }
    else if ( ! m_flattened_bulk_changes.empty() ) {
      encode_into(destination, m_flattened_bulk_changes);
    }
    else if ( ! m_nested_bulk_changes.empty() ) {
      encode_into(destination, m_nested_bulk_changes);
    }
    else if ( has_per_sample_changes() ) {
      encode_into(destination, m_per_sample_changes);
    }
    else {
      throw std::logic_error(std::string{"Encoding requested, but nothing to encode for individual field: "} + std::to_string(m_field_index));
    }
  }

 private:
  uint32_t m_field_index;                                                ///< index of this field in the Variant header
  int32_t m_field_type;                                                  ///< htslib type of this field, one of BCF_HT_*
  uint32_t m_num_samples;                                                ///< number of samples in this field
  uint32_t m_max_sample_value_length;                                    ///< length of the longest value for any sample in this field
  ENCODED_TYPE m_missing_value;                                          ///< BCF missing value for this field @note missing values are NOT constant per type. Eg., GT has different values from other int fields
  ENCODED_TYPE m_end_of_vector_value;                                    ///< BCF end-of-vector value for this field
  std::vector<BULK_CHANGE_TYPE> m_flattened_bulk_changes;                ///< stores bulk changes to this field in the form of a flattened and pre-padded vector
  std::vector<std::vector<BULK_CHANGE_TYPE>> m_nested_bulk_changes;      ///< stores bulk changes to this field in the form of a nested vector, with each inner vector representing values for one sample
  utils::ShortValueOptimizedStorage<ENCODED_TYPE> m_per_sample_changes;  ///< stores per-sample changes to this field efficiently
  bool m_removed;                                                        ///< keeps track of whether the caller has explicitly requested that this field be removed

  /**
   * Encoder functions for each type + storage combination.
   *
   * NOTE: non-virtual member functions are instantiated only if they are called,
   * so, eg., for the VariantBuilderIndivField<int32_t, int32_t> only the int32_t encode() methods
   * below will actually get instantiated. This turns out to be a pretty nice alternative to polymorphism.
   */

  void encode_into(kstring_t* destination, const std::vector<int32_t>& values) const {
    bcf_enc_int1(destination, m_field_index);

    // Since the user has provided us with a pre-flattened/padded set of values, we can just pass
    // it to htslib directly for encoding.
    // N.B. VariantBuilder has already ensured that values.size() is divisible by num_samples (unless validation is turned off)
    bcf_enc_vint(destination, values.size(), const_cast<int32_t*>(&(values[0])), values.size() / m_num_samples);
  }

  void encode_into(kstring_t* destination, const std::vector<std::vector<int32_t>>& values) const {
    auto min_value = INT32_MAX;
    auto max_value = INT32_MIN + 1;

    // We could just create a flattened vector and then pass it directly to htslib to spare us all this
    // low-level work, but it would be less efficient

    // First find the min/max values
    for ( const auto& sample_values : values ) {
      if ( sample_values.size() > 0 ) find_min_max_int_values(&(sample_values[0]), sample_values.size(), min_value, max_value);
    }

    // Then determine the resulting encoded type (int8, int16, or int32)
    const auto encoded_type = utils::int_encoded_type(min_value, max_value);

    // Encode the field index and type descriptor
    bcf_enc_int1(destination, m_field_index);
    bcf_enc_size(destination, m_max_sample_value_length, encoded_type);

    // Encode values for each sample, padding as necessary to the max length
    for ( const auto& sample_values : values ) {
      // OK to pass in nullptr to the encoder if a sample has no values
      encode_and_pad_int_values_as_type(destination, sample_values.empty() ? nullptr : &(sample_values[0]), sample_values.size(), m_max_sample_value_length, encoded_type);
    }
  }

  void encode_into(kstring_t* destination, const utils::ShortValueOptimizedStorage<int32_t>& values) const {
    auto max_length = values.max_value_length();
    auto min_value = INT32_MAX;
    auto max_value = INT32_MIN + 1;

    // We could just create a flattened vector and then pass it directly to htslib to spare us all this
    // low-level work, but it would be less efficient

    // First find the min/max values
    for ( auto i = 0u; i < values.capacity(); ++i ) {
      auto sample_value = values.get(i);
      if ( sample_value.second > 0 ) find_min_max_int_values(sample_value.first, sample_value.second, min_value, max_value);
    }

    // Then determine the resulting encoded type (int8, int16, or int32)
    const auto encoded_type = utils::int_encoded_type(min_value, max_value);

    // Encode field index and type descriptor
    bcf_enc_int1(destination, m_field_index);
    bcf_enc_size(destination, max_length, encoded_type);

    // Encode values for each sample, padding as necessary to the max length
    for ( auto i = 0u; i < values.capacity(); ++i ) {
      auto sample_value = values.get(i);
      // OK to pass in nullptr to the encoder if a sample has no values
      encode_and_pad_int_values_as_type(destination, sample_value.second == 0 ? nullptr : sample_value.first, sample_value.second, max_length, encoded_type);
    }
  }

  void encode_into(kstring_t* destination, const std::vector<float>& values) const {
    bcf_enc_int1(destination, m_field_index);

    // Since the user has provided us with a pre-flattened/padded set of values, we can just pass
    // it to htslib directly for encoding.
    // N.B. VariantBuilder has already ensured that values.size() is divisible by num_samples (unless validation is turned off)
    bcf_enc_size(destination, values.size() / m_num_samples, BCF_BT_FLOAT);
    kputsn(reinterpret_cast<const char*>(&(values[0])), values.size() * sizeof(float), destination);
  }

  void encode_into(kstring_t* destination, const std::vector<std::vector<float>>& values) const {
    // Encode field index and type descriptor
    bcf_enc_int1(destination, m_field_index);
    bcf_enc_size(destination, m_max_sample_value_length, BCF_BT_FLOAT);

    // Then encode values for each sample, padding as necessary to the max length
    for ( auto& sample_values : values ) {
      encode_and_pad_sample_values(destination, &(sample_values[0]), sample_values.size(), m_max_sample_value_length);
    }
  }

  void encode_into(kstring_t* destination, const utils::ShortValueOptimizedStorage<float>& values) const {
    const auto max_length = values.max_value_length();

    // Encode field index and type descriptor
    bcf_enc_int1(destination, m_field_index);
    bcf_enc_size(destination, max_length, BCF_BT_FLOAT);

    // Encode per-sample values, padding as necessary to the max length
    for ( auto i = 0u; i < values.capacity(); ++i ) {
      const auto sample_values = values.get(i);
      encode_and_pad_sample_values(destination, sample_values.first, sample_values.second, max_length);
    }
  }

  void encode_into(kstring_t* destination, const std::vector<std::string>& values) const {
    // Encode field index and type descriptor
    bcf_enc_int1(destination, m_field_index);
    bcf_enc_size(destination, m_max_sample_value_length, BCF_BT_CHAR);

    // Encode per-sample values, padding as necessary to the max length
    for ( const auto& str : values ) {
      encode_and_pad_sample_values(destination, str.c_str(), str.length(), m_max_sample_value_length);
    }
  }

  void encode_into(kstring_t* destination, const std::vector<std::vector<std::string>>& values) const {
    // this one is not possible, but is needed in order to allow the string instantiation of the template to compile
    throw std::logic_error("nested vectors of strings not supported");
  }

  void encode_into(kstring_t* destination, const utils::ShortValueOptimizedStorage<char>& values) const {
    const auto max_length = values.max_value_length();

    // Encode field index and type descriptor
    bcf_enc_int1(destination, m_field_index);
    bcf_enc_size(destination, max_length, BCF_BT_CHAR);

    // Encode per-sample values, padding as necessary to the max length
    for ( auto i = 0u; i < values.capacity(); ++i ) {
      const auto sample_values = values.get(i);
      encode_and_pad_sample_values(destination, sample_values.first, sample_values.second, max_length);
    }
  }

  /**
   * @brief Handles encoding and padding of per-sample values, with the exception of int8 and int16
   *        which are handled separately
   */
  void encode_and_pad_sample_values(kstring_t* destination, const ENCODED_TYPE* sample_values, const uint32_t num_values, const uint32_t field_width) const {
    // This is how much padding we usually need (exception is the missing sample case handled below)
    auto pad_size = field_width - num_values;

    // If this sample has values, copy them first
    if ( num_values > 0 ) {
      kputsn(reinterpret_cast<const char*>(sample_values), num_values * sizeof(ENCODED_TYPE), destination);
    }
    else {
      // If there are no values for this sample, we need to write a single missing value,
      // then pad with field_width - 1 end of vector values
      kputsn(reinterpret_cast<const char*>(&m_missing_value), sizeof(ENCODED_TYPE), destination);
      pad_size = field_width - 1;
    }

    // Pad with appropriate number of end of vector values
    pad_field(destination, reinterpret_cast<const char*>(&m_end_of_vector_value), sizeof(ENCODED_TYPE), pad_size);
  }

  /**
   * @brief Handles encoding and padding of integer per-sample values, branching based on the target type (int8, int16, or int32)
   */
  void encode_and_pad_int_values_as_type(kstring_t* destination, const int32_t* values, const uint32_t num_values, const uint32_t field_width, const uint8_t target_type) const {
    switch ( target_type ) {
      case BCF_BT_INT8:
        transcode_to_int8(destination, values, num_values, field_width);
        break;
      case BCF_BT_INT16:
        transcode_to_int16(destination, values, num_values, field_width);
        break;
      case BCF_BT_INT32:
        encode_and_pad_sample_values(destination, values, num_values, field_width);
        break;
      default:
        throw std::logic_error("Invalid target type in encode_and_pad_int_values_as_type()");
    }
  }

  /**
   * @brief A variation on the standard encode_and_pad_sample_values() function above that converts from int32
   *        to int8 as it encodes and pads the given sample values
   */
  void transcode_to_int8(kstring_t* destination, const int32_t* values, const uint32_t num_values, const uint32_t field_width) const {
    // missing value will vary depending on whether we're in the GT field or not
    const int8_t int8_missing = m_missing_value == bcf_int32_missing ? bcf_int8_missing : int8_t(m_missing_value);
    const int8_t int8_vector_end = bcf_int8_vector_end;  // need this to be an lvalue rather than a macro so we can take its address
    auto pad_size = field_width - num_values;

    // If there are values for this sample, convert them one at a time to int8
    if ( num_values > 0 ) {
      for ( auto i = 0u; i < num_values; ++i ) {
        if ( values[i] == m_end_of_vector_value ) kputc(int8_vector_end, destination);
        else if ( values[i] == m_missing_value ) kputc(int8_missing, destination);
        else kputc( values[i], destination);
      }
    }
    else {
      // If there are no values for this sample, we'll write an int8_missing followed by field_width - 1 int8_vector_ends
      kputc(int8_missing, destination);
      pad_size = field_width - 1;
    }

    // Pad field with appropriate number of vector end values
    pad_field(destination, reinterpret_cast<const char*>(&int8_vector_end), 1, pad_size);
  }

  /**
   * @brief A variation on the standard encode_and_pad_sample_values() function above that converts from int32
   *        to int16 as it encodes and pads the given sample values
   */
  void transcode_to_int16(kstring_t* destination, const int32_t* values, const uint32_t num_values, const uint32_t field_width) const {
    // missing value will vary depending on whether we're in the GT field or not
    const int16_t int16_missing = m_missing_value == bcf_int32_missing ? bcf_int16_missing : int16_t(m_missing_value);
    const int16_t int16_vector_end = bcf_int16_vector_end;  // need this to be an lvalue rather than a macro so we can take its address
    auto pad_size = field_width - num_values;
    int16_t int16_val = 0;

    // If there are values for this sample, convert them one at a time to int16
    if ( num_values > 0 ) {
      for ( auto i = 0u; i < num_values; ++i ) {
        if ( values[i] == m_end_of_vector_value ) int16_val = int16_vector_end;
        else if ( values[i] == m_missing_value ) int16_val = int16_missing;
        else int16_val = int16_t(values[i]);
        kputsn(reinterpret_cast<const char*>(&int16_val), 2, destination);
      }
    }
    else {
      // If there are no values for this sample, we'll write an int16_missing followed by field_width - 1 int16_vector_ends
      kputsn(reinterpret_cast<const char*>(&int16_missing), 2, destination);
      pad_size = field_width - 1;
    }

    // Pad field with appropriate number of vector end values
    pad_field(destination, reinterpret_cast<const char*>(&int16_vector_end), 2, pad_size);
  }

  void pad_field(kstring_t* destination, const char* padding, const uint32_t padding_bytes, const uint32_t pad_size) const {
    for ( auto i = 0u; i < pad_size; ++i ) {
      kputsn(padding, padding_bytes, destination);
    }
  }

  /**
   * @brief Locates the smallest and largest ints in values (ignoring missing/vector end values), and stores them
   *        in the provided min and max output parameters.
   */
  void find_min_max_int_values(const int32_t* values, const uint32_t num_values, int32_t& min, int32_t& max) const {
    // N.B. can't use std::minmax_element() here since we need to explicitly exclude the missing and vector end values
    std::for_each(values, values + num_values, [&min, &max] (int32_t val) {
      if ( val != bcf_int32_missing && val != bcf_int32_vector_end ) {
        if ( val > max ) max = val;
        if ( val < min ) min = val;
      }
    });
  }

  /**
   * @brief The length of the longest sample value in a flat, non-string vector is simply the length divided
   *        by the number of samples
   */
  template<class T>
  uint32_t max_sample_value_length(const std::vector<T>& sample_data) const {
    return sample_data.size() / m_num_samples;
  }

  /**
   * @brief The length of the longest sample value in a nested vector has to be determined manually
   */
  uint32_t max_sample_value_length(const std::vector<std::string>& sample_data) const {
    auto max_length = 0u;
    std::for_each(sample_data.begin(), sample_data.end(), [&max_length] (const std::string& str) {
      if ( str.length() > max_length ) max_length = str.length();
    });
    return max_length;
  }

  /**
   * @brief The length of the longest sample value in a flat vector of strings has to be determined manually
   */
  template<class T>
  uint32_t max_sample_value_length(const std::vector<std::vector<T>>& sample_data) const {
    auto max_length = 0u;
    std::for_each(sample_data.begin(), sample_data.end(), [&max_length] (const std::vector<T>& sample_values) {
      if ( sample_values.size() > max_length ) max_length = sample_values.size();
    });
    return max_length;
  }
};

}

#endif  /* gamgee__variant_builder_individual_field__guard */
