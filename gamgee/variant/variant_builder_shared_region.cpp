#include "variant_builder_shared_region.h"

#include "../missing.h"

#include <algorithm>
#include <stdexcept>
#include <stdlib.h>

using namespace std;

namespace gamgee {

VariantBuilderSharedRegion::VariantBuilderSharedRegion(const VariantHeader& header, const bool enable_validation) :
  m_header { header.m_header },   // Important: don't deep copy header here, call constructor that shares ownership
  m_shared_buffer { 0, 0, 0 },
  m_field_buffer_spans(uint32_t(SharedFieldIndex::INFO_START_INDEX) + header.n_shared_fields()),   // N.B.: calling constructor that takes an int, so use () instead of {}
  m_info_field_lookup_table(uint32_t(header.m_header->n[BCF_DT_ID])),   // N.B.: calling constructor that takes an int, so use () instead of {}
  m_ref_allele_length { 0 },
  m_num_alt_alleles { 0 },
  m_num_present_info_fields { 0 },
  m_enable_validation { enable_validation }
{
  // Our shared memory pool is allocated using malloc() so that htslib can freely realloc() it as we encode
  m_shared_buffer = utils::initialize_htslib_buffer(initial_shared_buffer_size);
  build_lookup_table();
  mark_all_fields_as_unset();
}

VariantBuilderSharedRegion::~VariantBuilderSharedRegion() {
  // Need to release the htslib malloc()'d memory we've been using as a buffer
  // (note: can't use smart pointers to manage the release of this memory for several good reasons,
  //  including the essential need to allow htslib to realloc() it freely)
  free(m_shared_buffer.s);
}

/**
 * @brief Create a table mapping from info field logical index (as declared in the VariantHeader) to physical
 *        offset from SharedFieldIndex::INFO_START_INDEX, which determines where the info field is located
 *        in m_field_buffer_spans
 */
void VariantBuilderSharedRegion::build_lookup_table() {
  auto info_field_count = 0;

  // Create an entry for all possible BCF_DT_ID indices in the VariantHeader
  for ( auto i = 0; i < m_header.m_header->n[BCF_DT_ID]; ++i ) {
    if ( m_header.has_shared_field(i) ) {
      m_info_field_lookup_table[i] = info_field_count++;
    }
    else {
      m_info_field_lookup_table[i] = unset_field_sentinel_value;
    }
  }
}

/**
 * @brief Mark a region of the shared buffer as in-use for storage of the encoded data for a particular field
 */
void VariantBuilderSharedRegion::mark_field_as_set(const uint32_t field_physical_index, const uint32_t previous_buffer_length) {
  // New field was encoded starting at offset previous_buffer_length in the shared buffer
  m_field_buffer_spans[field_physical_index].start_offset = previous_buffer_length;

  // Size of new field is the difference between the current buffer length and previous buffer length
  m_field_buffer_spans[field_physical_index].length = m_shared_buffer.l - previous_buffer_length;

  // If we've exceeded the threshold for acceptable buffer size, compact the buffer to eliminate
  // gaps caused by old/garbage data for fields
  if ( m_shared_buffer.l > buffer_compaction_threshold ) {
    compact_shared_buffer();
  }
}

void VariantBuilderSharedRegion::mark_all_fields_as_unset() {
  for ( auto field_index = 0u; field_index < m_field_buffer_spans.size(); ++field_index ) {
    mark_field_as_unset(field_index);
  }
}

void VariantBuilderSharedRegion::set_id(const std::string& id) {
  if ( missing(id) ) {
    remove_id();
    return;
  }

  encode_non_info_field(uint32_t(SharedFieldIndex::ID_INDEX), id);
}

void VariantBuilderSharedRegion::set_ref_allele(const std::string& ref_allele) {
  if ( m_enable_validation && missing(ref_allele) ) {
    throw invalid_argument{"Cannot set the required ref allele field to a missing value"};
  }

  encode_non_info_field(uint32_t(SharedFieldIndex::REF_ALLELE_INDEX), ref_allele);
  m_ref_allele_length = ref_allele.length();
}

void VariantBuilderSharedRegion::set_alt_allele(const std::string& alt_allele) {
  if ( missing(alt_allele) ) {
    remove_alt_alleles();
    return;
  }

  encode_non_info_field(uint32_t(SharedFieldIndex::ALT_ALLELES_INDEX), alt_allele);
  m_num_alt_alleles = 1;
}

void VariantBuilderSharedRegion::set_alt_alleles(const std::vector<std::string>& alt_alleles) {
  if ( alt_alleles.empty() ) {
    remove_alt_alleles();
    return;
  }

  if ( m_enable_validation ) {
    for ( const auto& allele : alt_alleles ) {
      if ( missing(allele) ) {
        throw invalid_argument{"set_alt_alleles(): missing alt allele encountered in non-empty vector of alt alleles"};
      }
    }
  }

  encode_non_info_field(uint32_t(SharedFieldIndex::ALT_ALLELES_INDEX), alt_alleles);
  m_num_alt_alleles = alt_alleles.size();
}

void VariantBuilderSharedRegion::set_filters(const std::vector<std::string>& filters) {
  if ( filters.empty() ) {
    remove_filters();
    return;
  }

  auto filter_ids = vector<int32_t>(filters.size());
  for ( auto i = 0u; i < filters.size(); ++i ) {
    filter_ids[i] = m_header.field_index(filters[i]);
    if ( m_enable_validation && ! m_header.has_filter(filter_ids[i]) ) {
      throw invalid_argument(string{"No filter with id "} + filters[i] + " in builder's header");
    }
  }

  encode_non_info_field(uint32_t(SharedFieldIndex::FILTERS_INDEX), filter_ids);
}

void VariantBuilderSharedRegion::set_filters(const std::vector<int32_t>& filters) {
  if ( filters.empty() ) {
    remove_filters();
    return;
  }

  if ( m_enable_validation ) {
    for ( const auto filter_id : filters ) {
      if ( ! m_header.has_filter(filter_id) ) {
        throw invalid_argument(string{"No filter with index "} + to_string(filter_id) + " in builder's header");
      }
    }
  }

  encode_non_info_field(uint32_t(SharedFieldIndex::FILTERS_INDEX), filters);
}

/**
 * @brief Produce a slight overestimate of the number of bytes required to hold all encoded data
 *        for this shared region
 */
uint32_t VariantBuilderSharedRegion::estimate_total_size() const {
  auto total_size = 0u;
  for ( const auto& field : m_field_buffer_spans ) {
    total_size += field.length;   // Unset/removed fields will have length 0
  }

  return total_size + 64;  // Overestimate a bit to increase chances of avoiding a realloc() in htslib
}

/**
 * @brief Copy all encoded data for all shared fields into the provided byte buffer in the proper order
 *        and format for final insertion into a Variant object.
 */
void VariantBuilderSharedRegion::encode_into(kstring_t* destination) const {
  // id and filters need explicit encoded missing values if absent,
  // ref allele is always present,
  // and a blank alt alleles field has no encoding at all

  if ( field_is_set(uint32_t(SharedFieldIndex::ID_INDEX)) ) {
    copy_encoded_field_into(destination, uint32_t(SharedFieldIndex::ID_INDEX));
  }
  else {
    bcf_enc_size(destination, 0, BCF_BT_CHAR);
  }

  if ( m_enable_validation && ! field_is_set(uint32_t(SharedFieldIndex::REF_ALLELE_INDEX)) ) {
    throw logic_error{"Missing required ref allele"};
  }
  copy_encoded_field_into(destination, uint32_t(SharedFieldIndex::REF_ALLELE_INDEX));

  if ( field_is_set(uint32_t(SharedFieldIndex::ALT_ALLELES_INDEX)) ) {
    copy_encoded_field_into(destination, uint32_t(SharedFieldIndex::ALT_ALLELES_INDEX));
  }

  if ( field_is_set(uint32_t(SharedFieldIndex::FILTERS_INDEX)) ) {
    copy_encoded_field_into(destination, uint32_t(SharedFieldIndex::FILTERS_INDEX));
  }
  else {
    bcf_enc_vint(destination, 0, 0, -1);
  }

  for ( auto info_field_index = uint32_t(SharedFieldIndex::INFO_START_INDEX); info_field_index < m_field_buffer_spans.size(); ++info_field_index ) {
    if ( field_is_set(info_field_index) ) {
      copy_encoded_field_into(destination, info_field_index);
    }
  }
}

/**
 * @brief Reset the shared region to a pristine state with no field data
 */
void VariantBuilderSharedRegion::clear() {
  // Mark buffer as empty, but reuse existing allocated capacity
  m_shared_buffer.l = 0;

  mark_all_fields_as_unset();

  m_ref_allele_length = 0;
  m_num_alt_alleles = 0;
  m_num_present_info_fields = 0;
}

void VariantBuilderSharedRegion::validate_info_field(const int32_t field_index, const uint32_t provided_type, const uint32_t num_values) const {
  validate_info_field_existence(field_index);

  if ( m_header.shared_field_type(field_index) != provided_type ) {
    throw invalid_argument(string{"Type mismatch for shared/info field with index "} + to_string(field_index));
  }

  // TODO: validate the cardinality of the field -- might need to be done at build time
}

void VariantBuilderSharedRegion::validate_info_field_existence(const int32_t field_index) const {
  if ( ! m_header.has_shared_field(field_index) ) {
    throw invalid_argument(string{"No shared/info field with index "} + to_string(field_index) + " found in builder's header");
  }
}

void VariantBuilderSharedRegion::encode_non_info_field(const uint32_t field_physical_index, const std::string& value) {
  const auto previous_buffer_length = m_shared_buffer.l;
  bcf_enc_vchar(&m_shared_buffer, value.length(), value.c_str());
  mark_field_as_set(field_physical_index, previous_buffer_length);
}

void VariantBuilderSharedRegion::encode_non_info_field(const uint32_t field_physical_index, const std::vector<std::string>& values) {
  const auto previous_buffer_length = m_shared_buffer.l;
  for ( const auto& str : values ) {
    bcf_enc_vchar(&m_shared_buffer, str.length(), str.c_str());
  }
  mark_field_as_set(field_physical_index, previous_buffer_length);
}

void VariantBuilderSharedRegion::encode_non_info_field(const uint32_t field_physical_index, const std::vector<int32_t>& values) {
  const auto previous_buffer_length = m_shared_buffer.l;
  bcf_enc_vint(&m_shared_buffer, values.size(), const_cast<int32_t*>(&(values[0])), -1);
  mark_field_as_set(field_physical_index, previous_buffer_length);
}

void VariantBuilderSharedRegion::encode_info_field(const int32_t field_logical_index, const uint32_t field_physical_index, const int32_t value) {
  const auto previous_buffer_length = m_shared_buffer.l;
  bcf_enc_int1(&m_shared_buffer, field_logical_index);
  bcf_enc_int1(&m_shared_buffer, value);
  mark_field_as_set(field_physical_index, previous_buffer_length);
}

void VariantBuilderSharedRegion::encode_info_field(const int32_t field_logical_index, const uint32_t field_physical_index, const std::vector<int32_t>& values) {
  const auto previous_buffer_length = m_shared_buffer.l;
  bcf_enc_int1(&m_shared_buffer, field_logical_index);
  bcf_enc_vint(&m_shared_buffer, values.size(), const_cast<int32_t*>(&(values[0])), -1);
  mark_field_as_set(field_physical_index, previous_buffer_length);
}

void VariantBuilderSharedRegion::encode_info_field(const int32_t field_logical_index, const uint32_t field_physical_index, const float value) {
  const auto previous_buffer_length = m_shared_buffer.l;
  bcf_enc_int1(&m_shared_buffer, field_logical_index);
  bcf_enc_vfloat(&m_shared_buffer, 1, const_cast<float*>(&value));
  mark_field_as_set(field_physical_index, previous_buffer_length);
}

void VariantBuilderSharedRegion::encode_info_field(const int32_t field_logical_index, const uint32_t field_physical_index, const std::vector<float>& values) {
  const auto previous_buffer_length = m_shared_buffer.l;
  bcf_enc_int1(&m_shared_buffer, field_logical_index);
  bcf_enc_vfloat(&m_shared_buffer, values.size(), const_cast<float*>(&(values[0])));
  mark_field_as_set(field_physical_index, previous_buffer_length);
}

void VariantBuilderSharedRegion::encode_info_field(const int32_t field_logical_index, const uint32_t field_physical_index, const std::string& value) {
  const auto previous_buffer_length = m_shared_buffer.l;
  bcf_enc_int1(&m_shared_buffer, field_logical_index);
  bcf_enc_vchar(&m_shared_buffer, value.length(), value.c_str());
  mark_field_as_set(field_physical_index, previous_buffer_length);
}

void VariantBuilderSharedRegion::encode_info_field(const int32_t field_logical_index, const uint32_t field_physical_index, const bool value) {
  const auto previous_buffer_length = m_shared_buffer.l;
  bcf_enc_int1(&m_shared_buffer, field_logical_index);
  bcf_enc_size(&m_shared_buffer, 0, BCF_BT_NULL);
  mark_field_as_set(field_physical_index, previous_buffer_length);
}

/**
 * @brief When the number of bytes in use in the shared buffer gets too large, we compact it to
 *        get rid of unused/garbage data and restore data locality
 */
void VariantBuilderSharedRegion::compact_shared_buffer() {
  // Don't compact the buffer if it's already compacted (this could theoretically happen if we have extremely
  // long field values such that the total size exceeds the buffer_compaction_threshold, but no space is actually
  // wasted in the buffer).
  if ( unused_buffer_space() == 0 ) {
    return;
  }

  // Otherwise, copy all encoded field values into a new buffer, with no gaps/wasted space between the values,
  // and update the m_field_buffer_spans mapping to match the new layout
  auto compacted_buffer = utils::initialize_htslib_buffer(initial_shared_buffer_size);

  for ( auto field_index = 0u; field_index < m_field_buffer_spans.size(); ++field_index ) {
    if ( field_is_set(field_index) ) {
      const auto prev_compacted_buffer_length = compacted_buffer.l;
      copy_encoded_field_into(&compacted_buffer, field_index);
      m_field_buffer_spans[field_index].start_offset = prev_compacted_buffer_length;
      m_field_buffer_spans[field_index].length = compacted_buffer.l - prev_compacted_buffer_length;
    }
  }

  // Release the old buffer and replace it with the newly compacted version
  free(m_shared_buffer.s);
  m_shared_buffer = compacted_buffer;
}

uint32_t VariantBuilderSharedRegion::unused_buffer_space() {
  auto used_space = 0u;
  for ( auto field_index = 0u; field_index < m_field_buffer_spans.size(); ++field_index ) {
    used_space += m_field_buffer_spans[field_index].length;
  }

  // Unused buffer space will be the different between the length of the current byte buffer (note: NOT
  // its capacity, which we don't care about when deciding whether to compact the buffer) minus the
  // buffer space actually used by all present fields.
  return m_shared_buffer.l - used_space;
}


}
