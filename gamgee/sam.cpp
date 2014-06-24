#include "htslib/sam.h"
#include "sam.h"
#include "utils/hts_memory.h"

#include <iostream>

using namespace std;

namespace gamgee {

/**
 * @brief creates a sam record that points to htslib memory already allocated
 *
 * @note the resulting Sam shares ownership of the pre-allocated memory via shared_ptr
 *       reference counting
 */
Sam::Sam(const std::shared_ptr<bam_hdr_t>& header, const std::shared_ptr<bam1_t>& body) noexcept :
  m_header {header},
  m_body {body}
{}

/**
 * @brief creates a deep copy of a sam record
 *
 * @note the copy will have exclusive ownership over the newly-allocated htslib memory
 *       until a data field (cigar, bases, etc.) is accessed, after which it will be
 *       shared via reference counting with the Cigar, etc. objects
 */
Sam::Sam(const Sam& other) :
  m_header { utils::make_shared_sam_header(utils::sam_header_deep_copy(other.m_header.get())) },
  m_body { utils::make_shared_sam(utils::sam_deep_copy(other.m_body.get())) }
{}

/**
 * @brief moves a sam record, transferring ownership of the underlying htslib memory
 */
Sam::Sam(Sam&& other) noexcept :
  m_header { move(other.m_header) },
  m_body { move(other.m_body) }
{}

/**
 * @brief creates a deep copy of a sam record
 *
 * @note the copy will have exclusive ownership over the newly-allocated htslib memory
 *       until a data field (cigar, bases, etc.) is accessed, after which it will be
 *       shared via reference counting with the Cigar, etc. objects
 */
Sam& Sam::operator=(const Sam& other) {
  if ( &other == this )  
    return *this;
  m_header = utils::make_shared_sam_header(utils::sam_header_deep_copy(other.m_header.get())); ///< shared_ptr assignment will take care of deallocating old sam record if necessary
  m_body = utils::make_shared_sam(utils::sam_deep_copy(other.m_body.get()));                   ///< shared_ptr assignment will take care of deallocating old sam record if necessary
  return *this;
}

/**
 * @brief moves a sam record, transferring ownership of the underlying htslib memory
 */
Sam& Sam::operator=(Sam&& other) noexcept {
  if ( &other == this ) 
    return *this;
  m_header = move(other.m_header); ///< shared_ptr assignment will take care of decrementing the reference count for the old managed object (and destroying it if necessary)
  m_body = move(other.m_body);     ///< shared_ptr assignment will take care of decrementing the reference count for the old managed object (and destroying it if necessary)
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
uint32_t Sam::unclipped_start() const {
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
uint32_t Sam::unclipped_stop() const {
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

/**
 * @brief retrieve an integer-valued tag by name
 *
 * @note returns a SamTag with is_present() == false if the read has no tag by this name
 */
SamTag<int32_t> Sam::integer_tag(const std::string& tag_name) const {
  const auto aux_ptr = bam_aux_get(m_body.get(), tag_name.c_str());
  return aux_ptr == nullptr ? SamTag<int32_t>(tag_name, 0, false) :
                              SamTag<int32_t>(tag_name, bam_aux2i(aux_ptr));

  // TODO: htslib bam_aux2i() returns 0 if the tag is not of integer type,
  //       but 0 is a valid value for a correctly-typed tag. This should be patched.
}

/**
 * @brief retrieve a double/float-valued tag by name
 *
 * @note returns a SamTag with is_present() == false if the read has no tag by this name
 */
SamTag<double> Sam::double_tag(const std::string& tag_name) const {
  const auto aux_ptr = bam_aux_get(m_body.get(), tag_name.c_str());
  return aux_ptr == nullptr ? SamTag<double>(tag_name, 0.0, false) :
                              SamTag<double>(tag_name, bam_aux2f(aux_ptr));

  // TODO: htslib bam_aux2f() returns 0.0 if the tag is not of double/float type,
  //       but 0.0 is a valid value for a correctly-typed tag. This should be patched.
}

/**
 * @brief retrieve a char-valued tag by name
 *
 * @note returns a SamTag with is_present() == false if the read has no tag by this name
 */
SamTag<char> Sam::char_tag(const std::string& tag_name) const {
  const auto aux_ptr = bam_aux_get(m_body.get(), tag_name.c_str());
  if ( aux_ptr == nullptr )  // tag doesn't exist
    return SamTag<char>(tag_name, '\0', false);

  const auto char_val = bam_aux2A(aux_ptr);
  if ( char_val == '\0' )    // tag not of char type
    return SamTag<char>(tag_name, '\0', false);

  return SamTag<char>(tag_name, char_val);
}

/**
 * @brief retrieve a string-valued tag by name
 *
 * @note returns a SamTag with is_present() == false if the read has no tag by this name
 */
SamTag<std::string> Sam::string_tag(const std::string& tag_name) const {
  const auto aux_ptr = bam_aux_get(m_body.get(), tag_name.c_str());
  if ( aux_ptr == nullptr )  // tag doesn't exist
    return SamTag<string>(tag_name, "", false);

  const auto str_ptr = bam_aux2Z(aux_ptr);
  if ( str_ptr == nullptr )  // tag not of string type
    return SamTag<string>(tag_name, "", false);

  return SamTag<string>(tag_name, string{str_ptr});
}

}

