#ifndef __gamgee__variant__
#define __gamgee__variant__

#include "variant_header.h"

#include "htslib/sam.h"

#include <string>
#include <memory>

namespace gamgee {

/**
 * @brief Utility class to manipulate a Variant record.
 */
class Variant {
 public:
  Variant() = default; ///< initializes a null header @warning if you need to create a Variant from scratch, use the builder instead
  explicit Variant(const std::shared_ptr<bcf_hdr_t>& header, const std::shared_ptr<bcf1_t>& body) noexcept;
  Variant(const Variant& other);
  Variant(Variant&& other) noexcept;
  Variant& operator=(const Variant& other);
  Variant& operator=(Variant&& other) noexcept;

  VariantHeader header() const { return VariantHeader{m_header}; }

  uint32_t chromosome()      const {return uint32_t(m_body->rid);}      ///< @brief returns the integer representation of the chromosome. Notice that chromosomes are listed in index order with regards to the header (so a 0-based number). Similar to Picards getReferenceIndex()
  uint32_t alignment_start() const {return uint32_t(m_body->pos+1);}    ///< @brief returns a 1-based alignment start position (as you would see in a VCF file). @note the internal encoding is 0-based to mimic that of the BCF files.
  uint32_t qual()            const {return uint32_t(m_body->qual);}     ///< @brief returns the Phred scaled site qual (probability that the site is not reference). See VCF spec.
  uint32_t n_samples()       const {return uint32_t(m_body->n_sample);} ///< @brief returns the number of samples in this Variant record
  uint32_t n_alleles()       const {return uint32_t(m_body->n_allele);} ///< @brief returns the number of alleles in this Variant record

  std::vector<uint8_t> genotype_quals() const; ///< @brief returns a vector with a copy of all the GQ values for all samples contiguously in memory.

 private:
  std::shared_ptr<bcf_hdr_t> m_header; ///< @brief htslib variant header pointer
  std::shared_ptr<bcf1_t> m_body;      ///< @brief htslib variant body pointer

  friend class VariantWriter;
};

}  // end of namespace

#endif /* defined(__gamgee__variant__) */
