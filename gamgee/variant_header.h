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
  VariantHeader() = default; ///< initializes a null header @warning if you need to create a VariantHeader from scratch, use the builder instead
  explicit VariantHeader(const std::shared_ptr<bcf_hdr_t>& header) : m_header{header} {} ///< @brief simple constructor that takes shared ownership of a header object
  VariantHeader(const VariantHeader& other);
  VariantHeader(VariantHeader&& other) noexcept;
  VariantHeader& operator=(const VariantHeader& other); 
  VariantHeader& operator=(VariantHeader&& other) noexcept; ///< @brief move assignment operator 
  ~VariantHeader() = default;

  /**
   * @brief returns the number of samples in the header 
   * @note much faster than getting the actual list of samples
   */
  uint32_t n_samples() const { return uint32_t(m_header->n[BCF_DT_SAMPLE]); }  

  std::vector<std::string> filters()       const; ///< @brief builds a vector with the filters 
  std::vector<std::string> samples()       const; ///< @brief builds a vector with the names of the samples 
  std::vector<std::string> contigs()       const; ///< @brief builds a vector with the contigs 
  std::vector<std::string> info_fields()   const; ///< @brief builds a vector with the info fields 
  std::vector<std::string> format_fields() const; ///< @brief builds a vector with the format fields

  void advanced_merge_header(const VariantHeader& other) { bcf_hdr_combine(other.m_header.get(), m_header.get()); }

 private:
  std::shared_ptr<bcf_hdr_t> m_header;
  
  friend class VariantWriter;
};

}

#endif
