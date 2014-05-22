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

  uint32_t chromosome()      const {return uint32_t(m_body->rid);}
  uint32_t alignment_start() const {return uint32_t(m_body->pos+1);}
  uint32_t qual()            const {return uint32_t(m_body->qual);}
  uint32_t n_samples()       const {return uint32_t(m_body->n_sample);}
  uint32_t n_alleles()       const {return uint32_t(m_body->n_allele);}

 private:
  std::shared_ptr<bcf_hdr_t> m_header;
  std::shared_ptr<bcf1_t> m_body;

  friend class VariantWriter;
};

}  // end of namespace

#endif /* defined(__gamgee__variant__) */
