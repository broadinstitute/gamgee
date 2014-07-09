#ifndef __gamgee__variant_header__
#define __gamgee__variant_header__

#include "htslib/vcf.h"

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
   * @brief returns the number of samples in the header 
   * @note much faster than getting the actual list of samples
   */
  uint32_t n_samples() const { return uint32_t(m_header->n[BCF_DT_SAMPLE]); }  

  std::vector<std::string> filters()       const; ///< @brief builds a vector with the filters 
  std::vector<std::string> samples()       const; ///< @brief builds a vector with the names of the samples 
  std::vector<std::string> chromosomes()   const; ///< @brief builds a vector with the contigs 
  std::vector<std::string> info_fields()   const; ///< @brief builds a vector with the info fields 
  std::vector<std::string> format_fields() const; ///< @brief builds a vector with the format fields
  bool has_format_field(const std::string)        const; ///< @brief checks if format field has the given field

  void advanced_merge_header(const VariantHeader& other) { bcf_hdr_combine(other.m_header.get(), m_header.get()); }

 private:
  std::shared_ptr<bcf_hdr_t> m_header;
  
  friend class VariantWriter;
};

}

#endif
