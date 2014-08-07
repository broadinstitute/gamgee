#ifndef gamgee__genotype_utils__guard
#define gamgee__genotype_utils__guard

#include "hts_memory.h"
#include "utils.h"
#include "format_field_type.h"
#include "../missing.h"

#include "htslib/vcf.h"

#include <memory>

namespace gamgee {

namespace utils {

using namespace std;

  /**
   * @brief Counts the genotype alleles.
   * @param format_ptr The GT field from the line.
   * @return the count of genotype alleles.
   * @warning format_ptr must be non-null.
   */
  inline uint32_t allele_count(const bcf_fmt_t* const format_ptr) {
    return format_ptr->n;
  }

  /**
   * @brief Returns true if the allele at position allele_index is missing.
   * @param format_ptr The GT field from the line.
   * @param data_ptr The GT for this sample.
   * @param allele_index The index within the GT field to return an integer allele.
   * @return true if the allele at position allele_index is missing.
   * @warning Only int8_t GT fields have been tested.
   * @warning Missing GT fields are untested.
   */
  bool allele_missing(const bcf_fmt_t* const format_ptr, const uint8_t* data_ptr, const uint32_t allele_index);

  template<class TYPE>
  inline int32_t allele_key(const uint8_t* data_ptr, const uint32_t allele_index, const TYPE missing) {
    // mostly copied from htslib
    const auto p = reinterpret_cast<const TYPE*>(data_ptr);
    if ( !(p[allele_index]>>1) || p[allele_index]==missing ) {
      return missing_values::int32;
    }
    return (p[allele_index]>>1)-1;
  }

  /**
   * @brief Returns the genotype allele at position allele_index.
   * @param format_ptr The GT field from the line.
   * @param data_ptr The GT for this sample.
   * @param allele_index The index within the GT field to return an integer allele.
   * @return The genotype allele at position allele_index.
   * @warning Only int8_t GT fields have been tested.
   */
  inline int32_t allele_key(
      const bcf_fmt_t* const format_ptr, const uint8_t* data_ptr, const uint32_t allele_index) {
    switch (format_ptr->type) {
    case BCF_BT_INT8:
      return allele_key<int8_t>(data_ptr, allele_index, bcf_int8_missing);
    case BCF_BT_INT16:
      return allele_key<int16_t>(data_ptr, allele_index, bcf_int16_missing);
    case BCF_BT_INT32:
      return allele_key<int32_t>(data_ptr, allele_index, bcf_int32_missing);
    default:
      throw invalid_argument("unknown GT field type: " + to_string(format_ptr->type));
    }
  }

  /**
   * @brief Returns the genotype allele keys.
   * @param body The shared memory variant "line" from a vcf, or bcf.
   * @param format_ptr The GT field from the line.
   * @param data_ptr The GT for this sample.
   * @return The genotype allele keys.
   * @warning Only int8_t GT fields have been tested.
   */
  vector<int32_t> allele_keys(const std::shared_ptr<bcf1_t>& body, const bcf_fmt_t* const format_ptr, const uint8_t* data_ptr);

  /**
   * @brief Returns the genotype allele strings.
   * @param body The shared memory variant "line" from a vcf, or bcf.
   * @param format_ptr The GT field from the line.
   * @param data_ptr The GT for this sample.
   * @return The genotype allele strings.
   * @warning Only int8_t GT fields have been tested.
   */
  vector<string> allele_strings(const std::shared_ptr<bcf1_t>& body,
      const bcf_fmt_t* const format_ptr, const uint8_t* data_ptr); // returns the actual alleles (e.g. A/T or TG/C or T/T/C or T/T/T/T/T ... )

  /**
   * @brief Returns the genotype allele string from this line.
   * @param body The shared memory variant "line" from a vcf, or bcf.
   * @param key_index The integer representation of the allele within this "line".
   * @return The genotype allele string from this line.
   * @warning Only int8_t GT fields have been tested.
   */
  string allele_key_to_string(const std::shared_ptr<bcf1_t>& body, const int32_t key_index);
}

}

#endif
