#ifndef __gamgee__sam__
#define __gamgee__sam__

#include "sam_header.h"
#include "sam_body.h"

#include "htslib/sam.h"

#include <string>

namespace gamgee {

/**
 * @brief Utility class to manipulate a Sam record.
 */
class Sam : public SamHeader , public SamBody {
 public:
  explicit Sam() = default;
  Sam(const bam_hdr_t* header, const bam1_t* body) : SamHeader{header} , SamBody{body} {}
};

}  // end of namespace

/**
 * @brief outputs the sam record in sam format.
 *
 * The output checks whether the record has quality scores. If it does, it outputs a sam record,
 * otherwise it outputs a fasta record.
 */
std::ostream& operator<< (std::ostream& os, const gamgee::Sam& sam);

#endif /* defined(__gamgee__sam__) */
