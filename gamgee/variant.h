#ifndef __gamgee__variant__
#define __gamgee__variant__

#include "variant_header.h"
#include "variant_body.h"

#include "htslib/sam.h"

#include <string>
#include <memory>

namespace gamgee {

/**
 * @brief Utility class to manipulate a Variant record.
 */
class Variant : public VariantBody {
 public:
  Variant() = default;
  explicit Variant(const std::shared_ptr<bcf1_t>& body, const std::shared_ptr<bcf_hdr_t>& header) noexcept : VariantBody{body}, m_header{header} {}
  VariantHeader header() {
    return VariantHeader{m_header}; ///< TODO: return a reference to a singleton VariantHeader object instead of constructing a new one each time, since we're already sharing the underlying htslib memory
  }

 private:
  std::shared_ptr<bcf_hdr_t> m_header;
};

}  // end of namespace

#endif /* defined(__gamgee__variant__) */
