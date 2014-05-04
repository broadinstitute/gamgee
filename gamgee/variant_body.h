#ifndef __gamgee__variant_body__
#define __gamgee__variant_body__

#include "htslib/vcf.h"

#include <string>
#include <memory>

namespace gamgee {

class VariantBody {
 public:
  VariantBody();
  explicit VariantBody(const std::shared_ptr<bcf1_t>& body);
  VariantBody(const VariantBody& other);
  VariantBody(VariantBody&& other) noexcept;
  VariantBody& operator=(const VariantBody& other);
  VariantBody& operator=(VariantBody&& other) noexcept;

  uint32_t chromosome()      const {return uint32_t(m_body->rid);}
  uint32_t alignment_start() const {return uint32_t(m_body->pos+1);}
  uint32_t qual()            const {return uint32_t(m_body->qual);}
  uint32_t n_samples()       const {return uint32_t(m_body->n_sample);}
  uint32_t n_alleles()       const {return uint32_t(m_body->n_allele);}

 private:
  std::shared_ptr<bcf1_t> m_body;

  friend class VariantWriter;
};

}
#endif
