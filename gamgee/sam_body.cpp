#include "htslib/sam.h"

#include "sam_body.h"

#include <iostream>

namespace gamgee {

/**
 * @brief creates an empty sam record and allocates memory for all the fields 
 * @note This constructor gives the object full ownership of the allocated memory. 
 */
SamBody::SamBody() :
  m_body {bam_init1()},
  m_must_destroy_body {true}
{}

/**
 * @brief creates a sam record that points to the htslib memory already allocated
 * @warning this constructor takes no ownership of the allocated data and is intended for higher performance in the usual case. If you need to take ownership of the data, call the copy_internal_record() member function.
 */
SamBody::SamBody(bam1_t* body) :
  m_body {body},
  m_must_destroy_body {false}
{}

/**
 * @brief deep copy of the given sam record
 * @note This constructor gives the object full ownership of the allocated memory. 
 */
SamBody::SamBody(const SamBody& other) :
  m_body {bam_init1()},
  m_must_destroy_body{true}
{
  bam_copy1(m_body, other.m_body);
}

/**
 * @brief takes over the ownership status of the other object as is (whether or not it owned anything)
 */
SamBody::SamBody(SamBody&& other) :
  m_body {other.m_body},
  m_must_destroy_body {other.m_must_destroy_body}
{
  other.m_body = nullptr;
}

/**
 * @brief deep copy of the given sam record
 * @note This constructor gives the object full ownership of the allocated memory and appropriately 
 * manages any already existing internal data.
 */
SamBody& SamBody::operator=(const SamBody& other) {
  if (!m_must_destroy_body)             ///< if this record was only pointing to someone else's body before, we need to allocate it's own
    copy_internal_record(other.m_body);
  else                                  ///< but if it already held it's own data, we can just copy over it
    bam_copy1(m_body, other.m_body);
  return *this;
}

/**
 * @brief destroys all htslib internal data if we own it
 */
SamBody::~SamBody() {
  if (m_must_destroy_body)
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

/** 
 * @brief takes ownership of the objects inside the record
 *
 * particularly useful if you are using an iterator (hence storing only weak
 * pointers to their data) and for some record you decide to keep it beyond
 * it's lifetime during the iteration. This function will create internal
 * copies of the data and destroy it accordingly when the object goes out of
 * scope (destructor)
 */
void SamBody::make_internal_copy() { 
  if (m_must_destroy_body)  
    return;                     ///< we already own it!
  copy_internal_record(m_body);
}

void SamBody::copy_internal_record(const bam1_t* record) {
  m_body = bam_init1();
  bam_copy1(m_body, record);
  m_must_destroy_body = true; ///< takes ownership of the memory and tells the destructor so
}


}

