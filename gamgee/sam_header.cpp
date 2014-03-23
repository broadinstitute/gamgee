#include "sam_header.h"

namespace gamgee {

  SamHeader::SamHeader() :
    m_header {bam_hdr_init()}
  {}

  SamHeader::SamHeader(const bam_hdr_t* header) :
    m_header {bam_hdr_dup(header)}
  {}

  SamHeader::SamHeader(const SamHeader& other) :
    m_header {bam_hdr_dup(other.m_header)}
  {}

  SamHeader::SamHeader(SamHeader&& other) :
    m_header {other.m_header}
  {
    other.m_header = nullptr;
  }

  SamHeader::~SamHeader() {
    bam_hdr_destroy(m_header);
  }

  SamHeader& SamHeader::operator=(const SamHeader& other) {
    // TODO: add check for self assignment!
    bam_hdr_destroy(m_header);
    m_header = bam_hdr_dup(other.m_header);
    return *this;
  }

}
