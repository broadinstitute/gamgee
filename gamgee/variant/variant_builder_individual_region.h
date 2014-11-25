#ifndef gamgee__variant_builder_individual_region__guard
#define gamgee__variant_builder_individual_region__guard

#include "variant.h"
#include "variant_builder_individual_field.h"

#include "htslib/kstring.h"
#include "htslib/vcf.h"

#include <vector>
#include <string>
#include <stdexcept>

namespace gamgee {

/**
 * @brief Helper class for VariantBuilder to manage the fields belonging to the individual region of Variant records.
 *
 * The individual region includes the various kinds of FORMAT fields.
 * This class manages the validation and bookkeeping, but not the encoding and storage, of these fields. The
 * storage/encoding for each field is handled by the lower-level VariantBuilderIndividualField class.
 */
class VariantBuilderIndividualRegion {
 public:
  explicit VariantBuilderIndividualRegion(const VariantHeader& header, const bool enable_validation);

  VariantBuilderIndividualRegion(VariantBuilderIndividualRegion&& other) = default;
  VariantBuilderIndividualRegion& operator=(VariantBuilderIndividualRegion&& other) = default;
  VariantBuilderIndividualRegion(const VariantBuilderIndividualRegion& other) = delete;
  VariantBuilderIndividualRegion& operator=(const VariantBuilderIndividualRegion& other) = delete;
  ~VariantBuilderIndividualRegion() = default;

  void set_enable_validation(const bool enable_validation) { m_enable_validation = enable_validation; }

  int32_t gt_index() const { return m_gt_field_index; }
  uint32_t num_present_fields() const { return m_num_present_fields; }
  bool modified() const { return m_num_present_fields > 0; }

  template<class FIELD_ID_TYPE, class BULK_FIELD_VALUES_TYPE>
  void bulk_set_genotype_field(const FIELD_ID_TYPE& field_id, BULK_FIELD_VALUES_TYPE&& field_values) {
    // Set boolean parameter in bulk_set_field() to true to indicate that we should allow GT to be set
    bulk_set_field(field_id, std::forward<BULK_FIELD_VALUES_TYPE>(field_values), BCF_HT_INT, m_int_fields, true);
  }

  template<class FIELD_ID_TYPE, class BULK_FIELD_VALUES_TYPE>
  void bulk_set_integer_field(const FIELD_ID_TYPE& field_id, BULK_FIELD_VALUES_TYPE&& field_values) {
    // Final boolean parameter to bulk_set_field() is kept as its default false here to disallow GT to be set as a regular int field
    bulk_set_field(field_id, std::forward<BULK_FIELD_VALUES_TYPE>(field_values), BCF_HT_INT, m_int_fields);
  }

  template<class FIELD_ID_TYPE, class BULK_FIELD_VALUES_TYPE>
  void bulk_set_float_field(const FIELD_ID_TYPE& field_id, BULK_FIELD_VALUES_TYPE&& field_values) {
    bulk_set_field(field_id, std::forward<BULK_FIELD_VALUES_TYPE>(field_values), BCF_HT_REAL, m_float_fields);
  }

  template<class FIELD_ID_TYPE, class BULK_FIELD_VALUES_TYPE>
  void bulk_set_string_field(const FIELD_ID_TYPE& field_id, BULK_FIELD_VALUES_TYPE&& field_values) {
    bulk_set_field(field_id, std::forward<BULK_FIELD_VALUES_TYPE>(field_values), BCF_HT_STR, m_string_fields);
  }

  template<class FIELD_ID_TYPE, class SAMPLE_ID_TYPE, class FIELD_VALUE_TYPE>
  void set_genotype_field_by_sample(const FIELD_ID_TYPE& field_id, const SAMPLE_ID_TYPE& sample_id, const FIELD_VALUE_TYPE* field_values, const uint32_t num_field_values) {
    // Set boolean parameter in set_field_by_sample() to true to indicate that we should allow GT to be set
    set_field_by_sample(field_id, sample_id, field_values, num_field_values, BCF_HT_INT, m_int_fields, true);
  }

  template<class FIELD_ID_TYPE, class SAMPLE_ID_TYPE, class FIELD_VALUE_TYPE>
  void set_integer_field_by_sample(const FIELD_ID_TYPE& field_id, const SAMPLE_ID_TYPE& sample_id, const FIELD_VALUE_TYPE* field_values, const uint32_t num_field_values) {
    // Final boolean parameter to set_field_by_sample() is kept as its default false here to disallow GT to be set as a regular int field
    set_field_by_sample(field_id, sample_id, field_values, num_field_values, BCF_HT_INT, m_int_fields);
  }

  template<class FIELD_ID_TYPE, class SAMPLE_ID_TYPE, class FIELD_VALUE_TYPE>
  void set_float_field_by_sample(const FIELD_ID_TYPE& field_id, const SAMPLE_ID_TYPE& sample_id, const FIELD_VALUE_TYPE* field_values, const uint32_t num_field_values) {
    set_field_by_sample(field_id, sample_id, field_values, num_field_values, BCF_HT_REAL, m_float_fields);
  }

  template<class FIELD_ID_TYPE, class SAMPLE_ID_TYPE, class FIELD_VALUE_TYPE>
  void set_string_field_by_sample(const FIELD_ID_TYPE& field_id, const SAMPLE_ID_TYPE& sample_id, const FIELD_VALUE_TYPE* field_values, const uint32_t num_field_values) {
    set_field_by_sample(field_id, sample_id, field_values, num_field_values, BCF_HT_STR, m_string_fields);
  }

  template<class FIELD_ID_TYPE>
  void remove_individual_field(const FIELD_ID_TYPE& field_id) {
    const auto field_idx = field_index(field_id);
    if ( m_enable_validation ) {
      validate_individual_field_existence(field_idx);
    }
    remove_field(field_idx);
  }

  uint32_t estimate_total_size() const;
  void encode_into(kstring_t* buffer) const;
  void clear();

 private:
  VariantHeader m_header;                                                         ///< header used for validation/field lookup purposes
  std::vector<int32_t> m_field_lookup_table;                                      ///< mapping from field logical index (as declared in the header) to physical location in the appropriate m_*_fields vector
  int32_t m_gt_field_index;                                                       ///< logical (header) index of the GT field
  uint32_t m_num_present_fields;                                                  ///< number of fields that have at least one non-missing value for a sample
  std::vector<VariantBuilderIndividualField<int32_t, int32_t>> m_int_fields;      ///< integer fields, indexed using m_field_lookup_table
  std::vector<VariantBuilderIndividualField<float, float>> m_float_fields;        ///< float fields, indexed using m_field_lookup_table
  std::vector<VariantBuilderIndividualField<char, std::string>> m_string_fields;  ///< string fields, indexed using m_field_lookup_table
  bool m_enable_validation;                                                       ///< should we validate?

  /**
   * Definitions of the size of a "short value" for each type -- used by the VariantBuilderIndividualField class
   * for storage optimization purposes.
   *
   * Note: can't use in-class static initialization here, as emplace_back() needs to take the address of these constants
   */
  static const uint32_t int_field_short_value_threshold;
  static const uint32_t float_field_short_value_threshold;
  static const uint32_t string_field_short_value_threshold;

  void build_lookup_tables();

  template<class FIELD_ID_TYPE, class BULK_FIELD_VALUES_TYPE, class FIELD_TYPE>
  void bulk_set_field(const FIELD_ID_TYPE& field_id, BULK_FIELD_VALUES_TYPE&& field_values, const int32_t provided_type, std::vector<FIELD_TYPE>& fields_of_type, const bool allow_gt = false) {
    const auto field_idx = field_index(field_id);
    if ( m_enable_validation ) {
      validate_individual_field(field_idx, provided_type, allow_gt);
      validate_multi_sample_vector_length(field_values);
    }

    auto& field = fields_of_type[m_field_lookup_table[field_idx]];
    const auto field_was_already_present = field.present();

    // Need to use std::forward() here so that we can handle both the move and the copy use cases
    field.set_entire_field(std::forward<BULK_FIELD_VALUES_TYPE>(field_values));

    // Field will not necessarily be present after setting (eg., field_values might have been something like { {}, {} }).
    // If the field is missing after giving it the user's value then treat it as an explicit user request to
    // remove the entire field.
    // Note that it's too expensive to check up front whether the value is something like { {}. {} }, which is why
    // we check afterwards instead.
    if ( ! field.present() ) {
      field.remove();
    }
    update_present_field_count(field_was_already_present, field.present());
  }

  template<class FIELD_ID_TYPE, class SAMPLE_ID_TYPE, class FIELD_VALUE_TYPE, class FIELD_TYPE>
  void set_field_by_sample(const FIELD_ID_TYPE& field_id, const SAMPLE_ID_TYPE& sample_id, const FIELD_VALUE_TYPE* field_values, const uint32_t num_field_values, const int32_t provided_type, std::vector<FIELD_TYPE>& fields_of_type, const bool allow_gt = false) {
    const auto field_idx = field_index(field_id);
    const auto sample_idx = sample_index(sample_id);
    if ( m_enable_validation ) {
      validate_individual_field(field_idx, sample_idx, provided_type, allow_gt);
    }

    auto& field = fields_of_type[m_field_lookup_table[field_idx]];
    const auto field_was_already_present = field.present();

    field.set_sample_field_value(sample_idx, field_values, num_field_values);

    // Field will not necessarily be present after setting (eg., field_values might have been empty).
    // Note that unlike in the bulk setting case, we don't treat empty values here as a user request to
    // remove the field, since we're only dealing with a single sample's data.
    update_present_field_count(field_was_already_present, field.present());
  }

  void remove_field(const int32_t field_index);

  // These overloads only exist to allow us to unify the string id / integer id cases in
  // the templated functions above
  int32_t field_index(const std::string& field_id) const { return m_header.field_index(field_id); }
  int32_t field_index(const uint32_t field_id) const { return int32_t(field_id); }
  int32_t sample_index(const std::string& sample_id) const { return m_header.sample_index(sample_id); }
  int32_t sample_index(const uint32_t sample_id) const { return int32_t(sample_id); }

  void validate_individual_field(const int32_t field_index, const uint32_t provided_type, const bool allow_gt) const;
  void validate_individual_field(const int32_t field_index, const int32_t sample_index, const uint32_t provided_type, const bool allow_gt) const;
  void validate_individual_field_existence(const int32_t field_index) const;

  template<class ELEMENT_TYPE>
  void validate_multi_sample_vector_length(const std::vector<std::vector<ELEMENT_TYPE>>& vec) const {
    // Empty vectors explicitly allowed; non-empty vectors must have a size EQUAL to the # of samples
    if ( vec.size() != m_header.n_samples() && ! vec.empty() ) {
      throw std::invalid_argument(std::string{"Number of elements in non-empty vector of vectors for individual field ("} + std::to_string(vec.size()) + ") not equal to the number of samples (" + std::to_string(m_header.n_samples()) + ")");
    }
  }

  template<class ELEMENT_TYPE>
  void validate_multi_sample_vector_length(const std::vector<ELEMENT_TYPE>& vec) const {
    const auto num_samples = m_header.n_samples();

    // Empty vectors explicitly allowed; non-empty vectors must have a size DIVISIBLE by # of samples
    if ( vec.size() % num_samples != 0 ) {
      throw std::invalid_argument(std::string{"Number of elements in flattened vector for individual field ("} + std::to_string(vec.size()) + ") not divisible by number of samples (" + std::to_string(num_samples) + ")");
    }
  }

  // need to specialize for the vector of string case
  void validate_multi_sample_vector_length(const std::vector<std::string>& vec) const {
    // Empty vectors explicitly allowed; non-empty vectors must have a size EQUAL to the # of samples
    if ( vec.size() != m_header.n_samples() && ! vec.empty() ) {
      throw std::invalid_argument(std::string{"Number of elements in non-empty vector for individual field ("} + std::to_string(vec.size()) + ") not equal to the number of samples (" + std::to_string(m_header.n_samples()) + ")");
    }
  }

  void update_present_field_count(const bool field_was_already_present, const bool field_currently_present) {
    if ( ! field_was_already_present && field_currently_present ) {
      ++m_num_present_fields;
    }
    else if ( field_was_already_present && ! field_currently_present ) {
      --m_num_present_fields;
    }
  }
};

}

#endif  /* gamgee__variant_builder_individual_region__guard */
