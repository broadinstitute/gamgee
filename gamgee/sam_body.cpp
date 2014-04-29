#include "htslib/sam.h"
#include "sam_body.h"
#include "hts_memory.h"

#include <iostream>

using namespace std;

namespace gamgee {

/**
 * @brief creates an empty sam record and allocates new htslib memory for all the fields
 *
 * @note the copy will have exclusive ownership over the newly-allocated htslib memory
 *       until a data field (cigar, bases, etc.) is accessed, after which it will be
 *       shared via reference counting with the Cigar, etc. objects
 */
SamBody::SamBody() :
  m_body {make_shared_bam(bam_init1())}
{}

/**
 * @brief creates a sam record that points to htslib memory already allocated
 *
 * @note the resulting SamBody shares ownership of the pre-allocated memory via shared_ptr
 *       reference counting
 */
SamBody::SamBody(const shared_ptr<bam1_t>& body) :
  m_body { body }
{}

/**
 * @brief creates a deep copy of a sam record
 *
 * @note the copy will have exclusive ownership over the newly-allocated htslib memory
 *       until a data field (cigar, bases, etc.) is accessed, after which it will be
 *       shared via reference counting with the Cigar, etc. objects
 */
SamBody::SamBody(const SamBody& other) :
  m_body { make_shared_bam(bam_deep_copy(other.m_body.get())) }
{}

/**
 * @brief moves a sam record, transferring ownership of the underlying htslib memory
 */
SamBody::SamBody(SamBody&& other) :
  m_body { move(other.m_body) }
{}

/**
 * @brief creates a deep copy of a sam record
 *
 * @note the copy will have exclusive ownership over the newly-allocated htslib memory
 *       until a data field (cigar, bases, etc.) is accessed, after which it will be
 *       shared via reference counting with the Cigar, etc. objects
 */
SamBody& SamBody::operator=(const SamBody& other) {
  if ( &other != this ) {
    // shared_ptr assignment will take care of decrementing the reference count for the
    // old managed object (and destroying it if necessary)
    m_body = make_shared_bam(bam_deep_copy(other.m_body.get()));
  }

  return *this;
}

/**
 * @brief moves a sam record, transferring ownership of the underlying htslib memory
 */
SamBody& SamBody::operator=(SamBody&& other) {
  if ( &other != this ) {
    // shared_ptr assignment will take care of decrementing the reference count for the
    // old managed object (and destroying it if necessary)
    m_body = move(other.m_body);
  }

  return *this;
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
  const auto* cigar = bam_get_cigar(m_body.get());
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
  const auto* cigar = bam_get_cigar(m_body.get());
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

