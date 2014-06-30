#ifndef gamgee__variant_builder_shared_region__guard
#define gamgee__variant_builder_shared_region__guard

#include "utils/hts_memory.h"
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "variant.h"

#include <string>
#include <vector>
#include <utility>

namespace gamgee {

/**
 * @brief Enum to represent the ordering of the various shared fields as they are physically laid out in the encoded data
 */
enum class SharedFieldIndex { ID_INDEX = 0, REF_ALLELE_INDEX = 1, ALT_ALLELES_INDEX = 2, FILTERS_INDEX = 3, INFO_START_INDEX = 4 };

/**
 * @brief Represents a section (range of bytes) in the shared memory pool VariantBuilderSharedRegion::m_shared_buffer
 *        that is currently in use to store encoded field data.
 */
struct SharedBufferSpan {
  int32_t start_offset;
  uint32_t length;
};

/**
 * @brief Helper class for VariantBuilder to manage the fields belonging to the shared region of Variant records.
 *
 * The shared region includes the ID, reference allele, alt alleles, filters, and INFO fields.
 * This class manages the validation, encoding, and storage of these fields.
 *
 * STORAGE STRATEGY: This class pre-allocates a large block of C memory (m_shared_buffer). When
 *                   we set a shared field we tell htslib to encode it directly into the end of
 *                   this buffer, keeping track of the buffer span used by each field. htslib
 *                   can resize the buffer if needed.
 *
 *                   Old field values remain in the buffer as unused sections (garbage). If too much
 *                   garbage accumulates we compact the buffer to restore locality, but in
 *                   general if the user of VariantBuilder is calling clear() regularly this
 *                   won't happen.
 *
 *                   When we encode the final shared byte buffer, we traverse the used regions of
 *                   m_shared_buffer to copy the encoded data for each field into the right place
 *                   in the final buffer.
 *
 *                   This approach to storage has the advantage of performing no allocations or
 *                   de-allocations in normal use, and most of the time the encoded data for all fields
 *                   will be adjacent in memory, albeit possibly out-of-order, allowing for very
 *                   cache-efficient encode operations.
 */
class VariantBuilderSharedRegion {
 public:
  explicit VariantBuilderSharedRegion(const VariantHeader& header, const bool enable_validation);

  VariantBuilderSharedRegion(VariantBuilderSharedRegion&& other) = default;
  VariantBuilderSharedRegion& operator=(VariantBuilderSharedRegion&& other) = default;
  VariantBuilderSharedRegion(const VariantBuilderSharedRegion& other) = delete;
  VariantBuilderSharedRegion& operator=(const VariantBuilderSharedRegion& other) = delete;
  ~VariantBuilderSharedRegion(); // N.B. non-default destruction!

  void set_enable_validation(const bool enable_validation) { m_enable_validation = enable_validation; }

  uint32_t ref_allele_length() const { return m_ref_allele_length; }
  uint32_t num_alt_alleles() const { return m_num_alt_alleles; }
  uint32_t num_present_info_fields() const { return m_num_present_info_fields; }
  bool modified() const {
    return field_is_set(uint32_t(SharedFieldIndex::ID_INDEX)) ||
           field_is_set(uint32_t(SharedFieldIndex::REF_ALLELE_INDEX)) ||
           field_is_set(uint32_t(SharedFieldIndex::ALT_ALLELES_INDEX)) ||
           field_is_set(uint32_t(SharedFieldIndex::FILTERS_INDEX)) ||
           m_num_present_info_fields > 0;
  }

  void set_id(const std::string& id);
  void set_ref_allele(const std::string& ref_allele);
  void set_alt_allele(const std::string& alt_allele);
  void set_alt_alleles(const std::vector<std::string>& alt_alleles);
  void set_filters(const std::vector<std::string>& filters);
  void set_filters(const std::vector<int32_t>& filters);

  template<class FIELD_ID_TYPE, class FIELD_VALUE_TYPE>
  void set_info_field(const FIELD_ID_TYPE& field_id, const FIELD_VALUE_TYPE& field_value, const int32_t provided_type, const uint32_t num_values) {
    const auto field_logical_index = info_field_logical_index(field_id);
    if ( m_enable_validation ) {
      validate_info_field(field_logical_index, provided_type, num_values);
    }
    const auto field_physical_index = info_field_physical_index(field_logical_index);
    const auto field_was_already_present = field_is_set(field_physical_index);

    // If the user passed in a missing value, treat it as a request to remove the field
    if ( missing(field_value) ) {
      mark_field_as_removed(field_physical_index);
      if ( field_was_already_present ) --m_num_present_info_fields;
    }
    else {
      encode_info_field(field_logical_index, field_physical_index, field_value);
      if ( ! field_was_already_present ) ++m_num_present_info_fields;
    }
  }

  void remove_id() {
    mark_field_as_removed(uint32_t(SharedFieldIndex::ID_INDEX));
  }

  void remove_alt_alleles() {
    mark_field_as_removed(uint32_t(SharedFieldIndex::ALT_ALLELES_INDEX));
    m_num_alt_alleles = 0;
  }

  void remove_filters() {
    mark_field_as_removed(uint32_t(SharedFieldIndex::FILTERS_INDEX));
  }

  template<class FIELD_ID_TYPE>
  void remove_info_field(const FIELD_ID_TYPE& field_id) {
    const auto field_logical_index = info_field_logical_index(field_id);
    if ( m_enable_validation ) {
      validate_info_field_existence(field_logical_index);
    }

    const auto field_physical_index = info_field_physical_index(field_logical_index);
    const auto field_was_already_present = field_is_set(field_physical_index);

    mark_field_as_removed(field_physical_index);
    if ( field_was_already_present ) --m_num_present_info_fields;
  }

  uint32_t estimate_total_size() const;
  void encode_into(kstring_t* destination) const;
  void clear();

 private:
  VariantHeader m_header;                                 ///< header used for validation/field lookup purposes
  kstring_t m_shared_buffer;                              ///< shared memory pool for storage of encoded data; allocated using malloc() so that htslib can realloc() it as needed
  std::vector<SharedBufferSpan> m_field_buffer_spans;     ///< mapping from field ids to used sections of the shared memory pool
  std::vector<uint32_t> m_info_field_lookup_table;        ///< mapping from field logical index to field physical index (ie., to the order in the m_field_buffer_spans)
  uint32_t m_ref_allele_length;                           ///< length of the reference allele
  uint32_t m_num_alt_alleles;                             ///< number of alt alleles
  uint32_t m_num_present_info_fields;                     ///< number of INFO fields set to a non-missing value
  bool m_enable_validation;                               ///< should we validate?

  static constexpr uint32_t initial_shared_buffer_size = 2000;    ///< initial capacity of m_shared_buffer (may increase over time)
  static constexpr uint32_t buffer_compaction_threshold = 10000;  ///< if usage of m_shared_buffer reaches this many bytes, we compact the memory pool to get rid of garbage/unused data
  static constexpr int32_t unset_field_sentinel_value = -1;       ///< value used to identify unset fields in m_field_buffer_spans
  static constexpr int32_t removed_field_sentinel_value = -2;     ///< value used to identify explicitly-removed fields in m_field_buffer_spans

  void build_lookup_table();

  bool field_is_set(const uint32_t field_physical_index) const {
    return ! field_is_unset(field_physical_index) && ! field_is_removed(field_physical_index);
  }

  void mark_field_as_set(const uint32_t field_physical_index, const uint32_t previous_buffer_length);

  bool field_is_unset(const uint32_t field_physical_index) const {
    return m_field_buffer_spans[field_physical_index].start_offset == unset_field_sentinel_value;
  }

  void mark_field_as_unset(const uint32_t field_physical_index) {
    m_field_buffer_spans[field_physical_index].start_offset = unset_field_sentinel_value;
    m_field_buffer_spans[field_physical_index].length = 0;
  }

  void mark_all_fields_as_unset();

  bool field_is_removed(const uint32_t field_physical_index) const {
    return m_field_buffer_spans[field_physical_index].start_offset == removed_field_sentinel_value;
  }

  void mark_field_as_removed(const uint32_t field_physical_index) {
    m_field_buffer_spans[field_physical_index].start_offset = removed_field_sentinel_value;
    m_field_buffer_spans[field_physical_index].length = 0;
  }

  uint32_t info_field_physical_index(const uint32_t field_logical_index) const {
    return uint32_t(SharedFieldIndex::INFO_START_INDEX) + m_info_field_lookup_table[field_logical_index];
  }

  // These trivial overloads only exist to allow us to unify the string id / integer id cases in
  // the templated info field setter functions above
  int32_t info_field_logical_index(const std::string& field_id) const { return m_header.field_index(field_id); }
  int32_t info_field_logical_index(const uint32_t field_id) const { return int32_t(field_id); }

  void validate_info_field(const int32_t field_index, const uint32_t provided_type, const uint32_t num_values) const;
  void validate_info_field_existence(const int32_t field_index) const;

  void encode_non_info_field(const uint32_t field_physical_index, const std::string& value);
  void encode_non_info_field(const uint32_t field_physical_index, const std::vector<std::string>& values);
  void encode_non_info_field(const uint32_t field_physical_index, const std::vector<int32_t>& values);

  void encode_info_field(const int32_t field_logical_index, const uint32_t field_physical_index, const int32_t value);
  void encode_info_field(const int32_t field_logical_index, const uint32_t field_physical_index, const std::vector<int32_t>& values);
  void encode_info_field(const int32_t field_logical_index, const uint32_t field_physical_index, const float value);
  void encode_info_field(const int32_t field_logical_index, const uint32_t field_physical_index, const std::vector<float>& values);
  void encode_info_field(const int32_t field_logical_index, const uint32_t field_physical_index, const std::string& value);
  void encode_info_field(const int32_t field_logical_index, const uint32_t field_physical_index, const bool value);

  void copy_encoded_field_into(kstring_t* target, const uint32_t field_physical_index) const {
    kputsn(m_shared_buffer.s + m_field_buffer_spans[field_physical_index].start_offset,
           int32_t(m_field_buffer_spans[field_physical_index].length),
           target);
  }

  void compact_shared_buffer();
  uint32_t unused_buffer_space();
};

}

#endif  /* gamgee__variant_builder_shared_region__guard */
