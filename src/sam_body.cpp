#include "htslib/sam.h"

#include "sam_body.h"

#include <iostream>

namespace gamgee {

  SamBody::SamBody() :
    m_body {bam_init1()}
  {}

  SamBody::SamBody(const bam1_t* body) :
    m_body {bam_init1()}
  {
    bam_copy1(m_body, body);
  }

  SamBody::SamBody(const SamBody& other) :
    m_body {bam_init1()}
  {
    bam_copy1(m_body, other.m_body);
  }

  SamBody::SamBody(SamBody&& other) :
    m_body {other.m_body}
  {
    other.m_body = nullptr;
  }

  SamBody& SamBody::operator=(const SamBody& other) {
    bam_copy1(m_body, other.m_body);
    return *this;
  }

  SamBody::~SamBody() {
    bam_destroy1(m_body);
  }

  /**
   * @brief calculates the theoretical alignment start of a read that has soft/hard-clips preceding the alignment
   *
   * @return the alignment start (1-based, inclusive) adjusted for clipped bases.  For example if the read
   * has an alignment start of 100 but the first 4 bases were clipped (hard or soft clipped)
   * then this method will return 96.
   *
   * Invalid to call on an unmapped read.
   * Invalid to call with cigar = null
   */
  uint32_t SamBody::unclipped_start() const {
    auto pos = alignment_start();
    const auto* cigar = bam_get_cigar(m_body);
    for (auto i = 0u; i != m_body->core.n_cigar; ++i) {
      const auto op = bam_cigar_op(cigar[i]);
      if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) 
        pos -= bam_cigar_oplen(cigar[i]);
      else
        break;
    }
    return pos;
  }

  /**
   * @brief calculates the theoretical alignment stop of a read that has soft/hard-clips preceding the alignment
   *
   * @return the alignment stop (1-based, inclusive) adjusted for clipped bases.  For example if the read
   * has an alignment stop of 100 but the last 4 bases were clipped (hard or soft clipped)
   * then this method will return 104.
   *
   * Invalid to call on an unmapped read.
   * Invalid to call with cigar = null
   */
  uint32_t SamBody::unclipped_stop() const {
    auto pos = alignment_stop();
    const auto* cigar = bam_get_cigar(m_body);
    for (auto i = int{m_body->core.n_cigar - 1}; i >= 0; --i) {
      const auto op = bam_cigar_op(cigar[i]);
      if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)
        pos += bam_cigar_oplen(cigar[i]);
      else
        break;
    }
    return pos;
  }


}
