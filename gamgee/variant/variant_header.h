#ifndef gamgee__variant_header__guard
#define gamgee__variant_header__guard

#include "htslib/vcf.h"

#include "../missing.h"

#include <memory>
#include <string>
#include <vector>

namespace gamgee {

template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
class VariantHeaderMerger;  //forward declaration to declare friendship in VariantHeader
/**
 * @brief Utility class to hold a variant header
 *
 * It can be used to read headers from a VCF/BCF file, but to create one from scratch 
 * you want to use the VariantHeaderBuilder
 *
 * Note: filters, shared fields, and individual fields all occupy the same index space, and users must be careful to access
 * these values appropriately.  As an example, consider one possible indexing scheme for a simple header with two of each:
 *
 * 0 FILTER_1
 * 1 SHARED_FIELD_1
 * 2 SHARED_FIELD_2
 * 3 INDIVIDUAL_FIELD_1
 * 4 FILTER_2
 * 5 INDIVIDUAL_FIELD_2
 *
 * Counters and string accessors work intuitively:
 *
 * n_filters() returns 2
 * n_shared_fields() returns 2
 * n_individual_fields() returns 2
 * filters() returns { "FILTER_1", "FILTER_2" }
 * shared_fields() returns { "SHARED_FIELD_1", "SHARED_FIELD_2" }
 * individual_fields() returns { "INDIVIDUAL_FIELD_1", "INDIVIDUAL_FIELD_2" }
 *
 * Iterating through these filters/fields requires more care.  While it's possible to use the string vectors for this,
 * use of those strings to retrieve fields will require an expensive field index lookup.  The counters are not usable
 * because they don't tell you which indices are for your desired type.
 *
 * Instead, use field_index_end() and the has_filter() / has_*_field() method appropriate for the type on each index.
 *
 * Sample usage:
 * for (auto idx = 0u; idx < header.field_index_end(); ++idx)
 *   if (header.has_individual_field(idx))
 *     do_something_with_field(variant.integer_individual_field(idx));
 */
class VariantHeader {
 public:
  VariantHeader() = default;                                                             ///< @brief initializes a null VariantHeader @warning if you need to create a VariantHeader from scratch, use the builder instead
  explicit VariantHeader(const std::shared_ptr<bcf_hdr_t>& header) : m_header{header} {} ///< @brief creates a VariantHeader given htslib object. @note used by all iterators
  VariantHeader(const VariantHeader& other);                                             ///< @brief makes a deep copy of a VariantHeader. Shared pointers maintain state to all other associated objects correctly.
  VariantHeader(VariantHeader&& other) noexcept;                                         ///< @brief moves VariantHeader accordingly. Shared pointers maintain state to all other associated objects correctly.
  VariantHeader& operator=(const VariantHeader& other);                                  ///< @brief deep copy assignment of a VariantHeader. Shared pointers maintain state to all other associated objects correctly.
  VariantHeader& operator=(VariantHeader&& other) noexcept;                              ///< @brief move assignment of a VariantHeader. Shared pointers maintain state to all other associated objects correctly.
  ~VariantHeader() = default;

  /**
   * @brief equality operators
   * @param rhs the other VariantHeader to compare to
   * @warn work in progress: does not fully compare underlying structures
   * @return whether these variant headers are equal
   */
  bool operator==(const VariantHeader& rhs) const;
  bool operator!=(const VariantHeader& rhs) const { return !operator==(rhs); }

  std::vector<std::string> samples() const;           ///< @brief builds a vector with the names of the samples
  uint32_t n_samples() const { return uint32_t(bcf_hdr_nsamples(m_header.get())); };  ///< @brief returns the number of samples in the header  @note much faster than getting the actual list of samples
  std::vector<std::string> chromosomes() const;       ///< @brief builds a vector with the contigs
  uint32_t n_chromosomes() const;                     ///< @brief returns the number of chromosomes declared in this header

  /**
  * @brief returns the last valid field index + 1, to indicate the end of field iteration
  *
  * Sample usage:
  * for (auto idx = 0u; idx < header.field_index_end(); ++idx)
  *   if (header.has_individual_field(idx))
  *     do_something_with_field(variant.integer_individual_field(idx));
  */
  uint32_t field_index_end() const { return static_cast<uint32_t>(m_header->n[BCF_DT_ID]); };

  /**
  * @brief returns the number of filters declared in this header
  * @warn do not use for iteration over filter indices -- use field_index_end() instead
  */
  uint32_t n_filters() const;
  std::vector<std::string> filters() const;           ///< @brief returns a vector of filter names

  /**
  * @brief returns the number of shared fields declared in this header
  * @warn do not use for iteration over filter indices -- use field_index_end() instead
  */
  uint32_t n_shared_fields() const;
  std::vector<std::string> shared_fields() const;     ///< @brief returns a vector of shared field names

  /**
  * @brief returns the number of individual fields declared in this header
  * @warn do not use for iteration over filter indices -- use field_index_end() instead
  */
  uint32_t n_individual_fields() const;
  std::vector<std::string> individual_fields() const; ///< @brief returns a vector of individual field names

  // type checking functions: returns BCF_HT_FLAG, BCF_HT_INT, BCF_HT_REAL, BCF_HT_STR from htslib/vcf.h
  uint8_t shared_field_type(const std::string& tag) const { return field_type(field_index(tag), BCF_HL_INFO); }     ///< @brief returns the type of this shared (INFO) field  @note must check whether the field exists before calling this function, as it doesn't check for you
  uint8_t shared_field_type(const int32_t index) const { return field_type(index, BCF_HL_INFO); }                   ///< @brief returns the type of this shared (INFO) field  @note must check whether the field exists before calling this function, as it doesn't check for you
  uint8_t individual_field_type(const std::string& tag) const { return field_type(field_index(tag), BCF_HL_FMT); }  ///< @brief returns the type of this individual (FORMAT) field  @note must check whether the field exists before calling this function, as it doesn't check for you
  uint8_t individual_field_type(const int32_t index) const { return field_type(index, BCF_HL_FMT); }                ///< @brief returns the type of this individual (FORMAT) field  @note must check whether the field exists before calling this function, as it doesn't check for you

  /**
   * returns the type of the field with the specified name and category (one of BCF_HL_FMT, BCF_HL_INFO, or BCF_HL_FLT)
   *
   * @note must check whether the field exists before calling this function, as it doesn't check for you
   */
  uint8_t field_type(const std::string& tag, const int32_t field_category) const { return field_type(field_index(tag), field_category); }

  /**
   * returns the type of the field with the specified index and category (one of BCF_HL_FMT, BCF_HL_INFO, or BCF_HL_FLT)
   *
   * @note must check whether the field exists before calling this function, as it doesn't check for you
   */
  uint8_t field_type(const int32_t index, const int32_t field_category) const { return bcf_hdr_id2type(m_header.get(), field_category, index); }
  /*
   * returns one of BCF_VL_* values for field with the specified name and category (one of BCF_HL_FMT, BCF_HL_INFO, or BCF_HL_FLT)
   *
   * @note must check whether the field exists before calling this function, as it doesn't check for you
   */
  uint32_t field_length_descriptor(const std::string& tag, const int32_t field_category) const { return field_length_descriptor(field_index(tag), field_category); }
  /**
   * returns one of BCF_VL_* values for the field with the specified index and category (one of BCF_HL_FMT, BCF_HL_INFO, or BCF_HL_FLT)
   *
   * @note must check whether the field exists before calling this function, as it doesn't check for you
   */
  uint32_t field_length_descriptor(const int32_t index, const int32_t field_category) const { return bcf_hdr_id2length(m_header.get(), field_category, index); }
  /**
   * returns number of values for the field with the specified name and category (one of BCF_HL_FMT, BCF_HL_INFO, or BCF_HL_FLT), 0xfffff for variable length fields
   *
   * @note must check whether the field exists before calling this function, as it doesn't check for you
   */
  uint32_t field_length(const std::string& tag, const int32_t field_category) const { return field_length(field_index(tag), field_category); }
  /**
   * returns number of values for the field with the specified index and category (one of BCF_HL_FMT, BCF_HL_INFO, or BCF_HL_FLT), 0xfffff for variable length fields
   *
   * @note must check whether the field exists before calling this function, as it doesn't check for you
   */
  uint32_t field_length(const int32_t index, const int32_t field_category) const { return bcf_hdr_id2number(m_header.get(), field_category, index); }
  /**
   * @brief checks whether the given filter is present given the filter name
   */
  bool has_filter(const std::string& filter_name) const { return has_field(field_index(filter_name), BCF_HL_FLT); }

  /**
   * @brief checks whether the given filter is present given the filter index
   */
  bool has_filter(const int32_t filter_index) const { return has_field(filter_index, BCF_HL_FLT); }

  /**
   * @brief checks whether the given shared (INFO) field is present given the field name
   */
  bool has_shared_field(const std::string& field_name) const { return has_field(field_index(field_name), BCF_HL_INFO); }

  /**
   * @brief checks whether the given shared (INFO) field is present given the field index
   */
  bool has_shared_field(const int32_t field_index) const { return has_field(field_index, BCF_HL_INFO); }

  /**
   * @brief checks whether the given individual (FORMAT) field is present given the field name
   */
  bool has_individual_field(const std::string& field_name) const { return has_field(field_index(field_name), BCF_HL_FMT); }

  /**
   * @brief checks whether the given individual (FORMAT) field is present given the field index
   */
  bool has_individual_field(const int32_t field_index) const { return has_field(field_index, BCF_HL_FMT); }

  /**
   * @brief checks whether the given field is present given the field name and field category (which must be one of BCF_HL_FMT, BCF_HL_INFO, or BCF_HL_FLT)
   */
  bool has_field(const std::string& field_name, const int32_t field_category) const { return has_field(field_index(field_name), field_category); }

  /**
   * @brief checks whether the given field is present given the field index and field category (one of BCF_HL_FMT, BCF_HL_INFO, or BCF_HL_FLT)
   */
  bool has_field(const int32_t field_index, const int32_t field_category) const {
    // Can't just use bcf_hdr_idinfo_exists() here since it assumes the index came from a hash lookup,
    // which is not always the case
    return field_index >= 0 &&
           field_index < m_header->n[BCF_DT_ID] &&
           m_header->id[BCF_DT_ID][field_index].val != nullptr &&
           m_header->id[BCF_DT_ID][field_index].val->hrec[field_category] != nullptr;
  }

  /**
   * @brief checks whether the given sample is present given the sample name
   */
  bool has_sample(const std::string& sample_name) const { return has_sample(sample_index(sample_name)); }

  /**
   * @brief checks whether the given sample is present given the sample index
   */
  bool has_sample(const int32_t sample_index) const {
    // Can't assume that sample_index came from a hash lookup, so must validate the hard way
    return sample_index >= 0 &&
           sample_index < m_header->n[BCF_DT_SAMPLE] &&
           m_header->id[BCF_DT_SAMPLE][sample_index].val != nullptr &&
           m_header->id[BCF_DT_SAMPLE][sample_index].val->id != -1;
  }

  /**
   * @brief looks up the index of a particular filter, shared or individual field tag, enabling subsequent O(1) random-access lookups for that field throughout the iteration. 
   * @return missing_values::int32_t if the tag is not present in the header (you can use missing() on the return value to check)
   * @note prefer this to looking up tag names during the iteration if you are looking for shared fields multiple times. 
   * @note if multiple fields (e.g. shared and individual) have the same tag (e.g. "DP"), they will also have the same index internally, so this function will do the right thing. The accessors for individual and shared field will know how to use the index to retrieve the correct field.
   */
  int32_t field_index(const std::string& tag) const {
    const auto index = bcf_hdr_id2int(m_header.get(), BCF_DT_ID, tag.c_str());
    return index >= 0 ? index : missing_values::int32;
  }

  /**
   * @brief looks up the index of a particular sample, enabling subsequent O(1) random-access lookups for that sample throughout the iteration.
   * @return missing_values::int32_t if the tag is not present in the header (you can use missing() on the return value to check)
   * @note prefer this to looking up sample names during the iteration if you are looking for samples multiple times.
   */
  int32_t sample_index(const std::string& sample) const {
    const auto index = bcf_hdr_id2int(m_header.get(), BCF_DT_SAMPLE, sample.c_str());
    return index >= 0 ? index : missing_values::int32;
  }

  std::string get_field_name(const int32_t field_idx) const {
    if(field_idx >= 0 && field_idx < m_header->n[BCF_DT_ID])
    {
      auto name_ptr = bcf_hdr_int2id(m_header.get(), BCF_DT_ID, field_idx);
      if(name_ptr)
	return name_ptr;
    }
    return "";
  }

  std::string get_sample_name(const int32_t sample_idx) const {
    if(sample_idx >= 0 && sample_idx < m_header->n[BCF_DT_SAMPLE])
    {
      auto name_ptr= bcf_hdr_int2id(m_header.get(), BCF_DT_SAMPLE, sample_idx);
      if(name_ptr)
	return name_ptr;
    }
    return "";
  }

  std::vector<bcf_hrec_t*> advanced_all_header_fields()
      const;  ///< Retrieve all the header fields with all its key/value
              ///parings. Lifetime of the pointers is the same as this header
              ///object. This is an advanced API option, if you are using it you
              ///should know what you are doing.

 private:
  std::shared_ptr<bcf_hdr_t> m_header;

  friend class Variant;
  friend class VariantWriter;
  friend class VariantHeaderBuilder;
  friend class VariantBuilder;       ///< builder needs access to the internals in order to build efficiently
  friend class VariantBuilderSharedRegion;
  friend class VariantBuilderIndividualRegion;
  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  friend class VariantHeaderMerger;  //to access m_header
};

}

#endif // gamgee__variant_header__guard
