#ifndef gamgee__variant_header__guard
#define gamgee__variant_header__guard

#include "htslib/vcf.h"
#include "missing.h"

#include <memory>
#include <string>
#include <vector>

namespace gamgee {
  
/**
 * @brief Utility class to hold a variant header
 *
 * It can be used to read headers from a VCF/BCF file, but to create one from scratch 
 * you want to use the VariantHeaderBuilder
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

  std::vector<std::string> filters() const;           ///< @brief builds a vector with the filters
  uint32_t n_filters() const;                         ///< @brief returns the number of filters declared in this header
  std::vector<std::string> samples() const;           ///< @brief builds a vector with the names of the samples
  uint32_t n_samples() const { return uint32_t(bcf_hdr_nsamples(m_header.get())); };  ///< @brief returns the number of samples in the header  @note much faster than getting the actual list of samples
  std::vector<std::string> chromosomes() const;       ///< @brief builds a vector with the contigs
  uint32_t n_chromosomes() const;                     ///< @brief returns the number of chromosomes declared in this header
  std::vector<std::string> shared_fields() const;     ///< @brief builds a vector with the info fields
  uint32_t n_shared_fields() const;                   ///< @brief returns the number of shared fields declared in this header
  std::vector<std::string> individual_fields() const; ///< @brief builds a vector with the format fields
  uint32_t n_individual_fields() const;               ///< @brief returns the number of individual fields declared in this header

  // type checking functions: returns BCF_HT_FLAG, BCF_HT_INT, BCF_HT_REAL, BCF_HT_STR from htslib/vcf.h
  uint8_t shared_field_type(const std::string& tag) const;      ///< returns the type of this shared (INFO) field
  uint8_t shared_field_type(const int32_t index) const;         ///< returns the type of this shared (INFO) field
  uint8_t individual_field_type(const std::string& tag) const;  ///< returns the type of this individual (FORMAT) field
  uint8_t individual_field_type(const int32_t index) const;     ///< returns the type of this individual (FORMAT) field

  /**
   * @brief checks whether the given filter is present given the filter name
   */
  bool has_filter(const std::string& filter_name) const { return field_exists(field_index(filter_name), BCF_HL_FLT); }

  /**
   * @brief checks whether the given filter is present given the filter index
   */
  bool has_filter(const int32_t filter_index) const { return field_exists(filter_index, BCF_HL_FLT); }

  /**
   * @brief checks whether the given shared (INFO) field is present given the field name
   */
  bool has_shared_field(const std::string& field_name) const { return field_exists(field_index(field_name), BCF_HL_INFO); }

  /**
   * @brief checks whether the given shared (INFO) field is present given the field index
   */
  bool has_shared_field(const int32_t field_index) const { return field_exists(field_index, BCF_HL_INFO); }

  /**
   * @brief checks whether the given individual (FORMAT) field is present given the field name
   */
  bool has_individual_field(const std::string& field_name) const { return field_exists(field_index(field_name), BCF_HL_FMT); }

  /**
   * @brief checks whether the given individual (FORMAT) field is present given the field index
   */
  bool has_individual_field(const int32_t field_index) const { return field_exists(field_index, BCF_HL_FMT); }

  /**
   * @brief checks whether the given sample is present given the sample name
   */
  bool has_sample(const std::string& sample_name) const { return sample_exists(sample_index(sample_name)); }

  /**
   * @brief checks whether the given sample is present given the sample index
   */
  bool has_sample(const int32_t sample_index) const { return sample_exists(sample_index); }

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

 private:
  std::shared_ptr<bcf_hdr_t> m_header;


  bool field_exists(const int32_t field_index, const int32_t field_type) const {
    // Can't just use bcf_hdr_idinfo_exists() here since it assumes the index came from a hash lookup,
    // which is not always the case
    return field_index >= 0 &&
           field_index < m_header->n[BCF_DT_ID] &&
           m_header->id[BCF_DT_ID][field_index].val != nullptr &&
           m_header->id[BCF_DT_ID][field_index].val->hrec[field_type] != nullptr;
  }

  bool sample_exists(const int32_t sample_index) const {
    // Can't assume that sample_index came from a hash lookup, so must validate the hard way
    return sample_index >= 0 &&
           sample_index < m_header->n[BCF_DT_SAMPLE] &&
           m_header->id[BCF_DT_SAMPLE][sample_index].val != nullptr &&
           m_header->id[BCF_DT_SAMPLE][sample_index].val->id != -1;
  }
  

  friend class VariantWriter;
  friend class VariantHeaderBuilder;
  friend class VariantBuilder;       ///< builder needs access to the internals in order to build efficiently
  friend class VariantBuilderSharedRegion;
  friend class VariantBuilderIndividualRegion;
};

}

#endif // gamgee__variant_header__guard
