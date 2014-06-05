#ifndef __gamgee__variant__
#define __gamgee__variant__

#include "variant_header.h"
#include "variant_field.h"
#include "variant_field_value.h"


#include "htslib/sam.h"

#include <string>
#include <memory>

namespace gamgee {

/**
 * @brief Utility class to manipulate a Variant record.
 */
class Variant {
 public:
  Variant() = default;                                                                                      ///< @brief initializes a null Variant @note this is only used internally by the iterators @warning if you need to create a Variant from scratch, use the builder instead
  explicit Variant(const std::shared_ptr<bcf_hdr_t>& header, const std::shared_ptr<bcf1_t>& body) noexcept; ///< @brief creates a Variant given htslib objects. @note used by all iterators
  Variant(const Variant& other);                                                                            ///< @brief makes a deep copy of a Variant and it's header. Shared pointers maintain state to all other associated objects correctly.
  Variant(Variant&& other) noexcept;                                                                        ///< @brief moves Variant and it's header accordingly. Shared pointers maintain state to all other associated objects correctly.
  Variant& operator=(const Variant& other);                                                                 ///< @brief deep copy assignment of a Variant and it's header. Shared pointers maintain state to all other associated objects correctly.
  Variant& operator=(Variant&& other) noexcept;                                                             ///< @brief move assignment of a Variant and it's header. Shared pointers maintain state to all other associated objects correctly.

  VariantHeader header() const { return VariantHeader{m_header}; }

  uint32_t chromosome()      const {return uint32_t(m_body->rid);}      ///< @brief returns the integer representation of the chromosome. Notice that chromosomes are listed in index order with regards to the header (so a 0-based number). Similar to Picards getReferenceIndex()
  uint32_t alignment_start() const {return uint32_t(m_body->pos+1);}    ///< @brief returns a 1-based alignment start position (as you would see in a VCF file). @note the internal encoding is 0-based to mimic that of the BCF files.
  uint32_t qual()            const {return uint32_t(m_body->qual);}     ///< @brief returns the Phred scaled site qual (probability that the site is not reference). See VCF spec.
  uint32_t n_samples()       const {return uint32_t(m_body->n_sample);} ///< @brief returns the number of samples in this Variant record
  uint32_t n_alleles()       const {return uint32_t(m_body->n_allele);} ///< @brief returns the number of alleles in this Variant record

  VariantField<VariantFieldValue<int32_t>> genotype_quals() const;   ///< @brief returns an iterable object with a pointer to all the GQ values for all samples contiguously in memory.
  VariantField<VariantFieldValue<int32_t>> phred_likelihoods() const;


 private:
  std::shared_ptr<bcf_hdr_t> m_header; ///< @brief htslib variant header pointer
  std::shared_ptr<bcf1_t> m_body;      ///< @brief htslib variant body pointer

  friend class VariantWriter;
};

}  // end of namespace

#endif /* defined(__gamgee__variant__) */
