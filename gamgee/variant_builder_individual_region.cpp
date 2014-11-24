#include "variant_builder_individual_region.h"

using namespace std;

namespace gamgee {

/**
 * Definitions of the size of a "short value" for each type -- used by the VariantBuilderIndividualField class
 * for storage optimization purposes.
 */
const uint32_t VariantBuilderIndividualRegion::int_field_short_value_threshold = 5;      // <= 5 ints = a "short" int field for storage purposes
const uint32_t VariantBuilderIndividualRegion::float_field_short_value_threshold = 5;    // <= 5 floats = a "short" float field for storage purposes
const uint32_t VariantBuilderIndividualRegion::string_field_short_value_threshold = 20;  // <= 20 chars = a "short" string field for storage purposes


VariantBuilderIndividualRegion::VariantBuilderIndividualRegion(const VariantHeader& header, const bool enable_validation):
    m_header { header.m_header },   // Important: don't deep copy header here, call constructor that shares ownership
    m_field_lookup_table(uint32_t(header.m_header->n[BCF_DT_ID])),   // N.B.: calling constructor that takes an int, so use () instead of {}
    m_gt_field_index { header.field_index("GT") },
    m_num_present_fields { 0 },
    m_int_fields{},
    m_float_fields{},
    m_string_fields{},
    m_enable_validation { enable_validation }
{
  const auto num_indiv_fields = header.n_individual_fields();
  m_int_fields.reserve(num_indiv_fields);
  m_float_fields.reserve(num_indiv_fields);
  m_string_fields.reserve(num_indiv_fields);
  build_lookup_tables();
}

/**
 * @brief Construct individual field objects and create a mapping from field logical (header) index to
 *        physical locations in the m_*_fields vectors.
 */
void VariantBuilderIndividualRegion::build_lookup_tables() {
  // Must copy float missing/eov values into local variables; using bcf_float_missing and bcf_float_vector_end
  // values directly doesn't work
  auto float_missing = 0.0f; bcf_float_set_missing(float_missing);
  auto float_vector_end = 0.0f; bcf_float_set_vector_end(float_vector_end);
  const auto num_samples = m_header.n_samples();

  // Create an entry for all possible BCF_DT_ID indices in the VariantHeader
  for ( auto i = 0; i < m_header.m_header->n[BCF_DT_ID]; ++i ) {
    if ( m_header.has_individual_field(i) ) {
      const auto field_type = m_header.individual_field_type(i);
      switch ( field_type ) {
        case BCF_HT_INT :
          m_int_fields.emplace_back(num_samples, i, field_type, bcf_int32_missing, bcf_int32_vector_end, int_field_short_value_threshold);
          m_field_lookup_table[i] = m_int_fields.size() - 1;
          break;
        case BCF_HT_REAL :
          m_float_fields.emplace_back(num_samples, i, field_type, float_missing, float_vector_end, float_field_short_value_threshold);
          m_field_lookup_table[i] = m_float_fields.size() - 1;
          break;
        case BCF_HT_STR :
          // GT is a string field in the header, although we encode it as an int field
          if ( i == m_gt_field_index ) {
            // Note that the missing value for genotypes is 0, not bcf_int32_missing
            m_int_fields.emplace_back(num_samples, i, field_type, 0, bcf_int32_vector_end, int_field_short_value_threshold);
            m_field_lookup_table[i] = m_int_fields.size() - 1;
          }
          else {
            // TODO: should we use bcf_str_missing here? htslib generally uses '.' instead of it...
            m_string_fields.emplace_back(num_samples, i, field_type, '.', bcf_str_vector_end, string_field_short_value_threshold);
            m_field_lookup_table[i] = m_string_fields.size() - 1;
          }
          break;
        default:
          throw logic_error(string{"Found format field in header with unsupported type: "} + to_string(field_type));
      }
    }
    else {
      m_field_lookup_table[i] = -1;
    }
  }
}

/**
 * @brief Produce a slight overestimate of the total size of the encoded data for this individual region
 */
uint32_t VariantBuilderIndividualRegion::estimate_total_size() const {
  auto buffer_size = 0u;

  for_each(m_int_fields.begin(), m_int_fields.end(), [&buffer_size, this] (const VariantBuilderIndividualField<int32_t, int32_t>& field) {
    buffer_size += field.estimated_encoded_size();
  });
  for_each(m_float_fields.begin(), m_float_fields.end(), [&buffer_size, this] (const VariantBuilderIndividualField<float, float>& field) {
    buffer_size += field.estimated_encoded_size();
  });
  for_each(m_string_fields.begin(), m_string_fields.end(), [&buffer_size, this] (const VariantBuilderIndividualField<char, string>& field) {
    buffer_size += field.estimated_encoded_size();
  });

  return buffer_size + 64;  // intentionally overestimate to try to avoid a realloc() later
}

/**
 * @brief Encode all individual fields into the provided byte buffer in the proper order
 *        and format for final insertion into a Variant object.
 */
void VariantBuilderIndividualRegion::encode_into(kstring_t* buffer) const {
  // Theoretically possible that GT is not declared at all in the header
  const auto gt_field_is_declared = ! missing(m_gt_field_index);
  const auto gt_physical_index = gt_field_is_declared ? m_field_lookup_table[m_gt_field_index] : missing_values::int32;

  // GT must be encoded first, as per the spec
  if ( gt_field_is_declared && m_int_fields[gt_physical_index].present() ) {
    m_int_fields[gt_physical_index].encode_into(buffer);
  }

  // TODO: order the remaining format fields in a more sensible (or at least customizable) way than by type,
  //       such as by global field index
  for ( auto i = 0u; i < m_int_fields.size(); ++i ) {
    if ( (int32_t(i) != gt_physical_index || ! gt_field_is_declared) && m_int_fields[i].present() ) {
      m_int_fields[i].encode_into(buffer);
    }
  }

  for ( auto& field : m_float_fields ) {
    if ( field.present() ) field.encode_into(buffer);
  }

  for ( auto& field : m_string_fields ) {
    if ( field.present() ) field.encode_into(buffer);
  }
}

/**
 * @brief Reset the individual region to a pristine state with no field data
 */
void VariantBuilderIndividualRegion::clear() {
  for ( auto& field : m_int_fields ) { field.clear(); }
  for ( auto& field : m_float_fields ) { field.clear(); }
  for ( auto& field : m_string_fields ) { field.clear(); }
  m_num_present_fields = 0;
}

void VariantBuilderIndividualRegion::remove_field(const int32_t field_index) {
  const auto field_type = field_index == m_gt_field_index ? BCF_HT_INT : m_header.individual_field_type(field_index);
  auto field_was_already_present = false;

  switch ( field_type ) {
    case BCF_HT_INT:
      field_was_already_present = m_int_fields[m_field_lookup_table[field_index]].present();
      m_int_fields[m_field_lookup_table[field_index]].remove();
      break;
    case BCF_HT_REAL:
      field_was_already_present = m_float_fields[m_field_lookup_table[field_index]].present();
      m_float_fields[m_field_lookup_table[field_index]].remove();
      break;
    case BCF_HT_STR:
      field_was_already_present = m_string_fields[m_field_lookup_table[field_index]].present();
      m_string_fields[m_field_lookup_table[field_index]].remove();
      break;
  }

  if ( field_was_already_present ) {
    --m_num_present_fields;
  }
}

void VariantBuilderIndividualRegion::validate_individual_field(const int32_t field_index, const uint32_t provided_type, const bool allow_gt) const {
  validate_individual_field_existence(field_index);

  // Special-case the GT field when doing type checking, since its nominal type in the header will be BCF_HT_STR,
  // but we require integer data encoded via one of the genotype setter functions
  if ( field_index == m_gt_field_index ) {
    if ( ! allow_gt || provided_type != BCF_HT_INT ) throw invalid_argument(string{"Type mismatch for GT field: must set GT using a genotype-specific setter function, and provide integer data"});
  }
  else if ( m_header.individual_field_type(field_index) != provided_type ) {
    throw invalid_argument(string{"Type mismatch for individual field with index "} + to_string(field_index));
  }

  // TODO: validate cardinality of the fields (must be done at build time, and could be expensive...)
}

void VariantBuilderIndividualRegion::validate_individual_field(const int32_t field_index, const int32_t sample_index, const uint32_t provided_type, const bool allow_gt) const {
  validate_individual_field(field_index, provided_type, allow_gt);

  if ( ! m_header.has_sample(sample_index) ) {
    throw invalid_argument(string{"No sample with index "} + to_string(sample_index) + " found in builder's header");
  }
}

void VariantBuilderIndividualRegion::validate_individual_field_existence(const int32_t field_index) const {
  if ( ! m_header.has_individual_field(field_index) ) {
    throw invalid_argument(string{"No individual field with index "} + to_string(field_index) + " found in builder's header");
  }
}

}
