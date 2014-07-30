#include "sam_header.h"
#include "utils/hts_memory.h"

using namespace std;

namespace gamgee {

  /**
   * @brief creates a SamHeader object that points to htslib memory already allocated
   *
   * @note the resulting SamHeader object shares ownership of the pre-allocated memory via
   *       shared_ptr reference counting
   */
  SamHeader::SamHeader(const std::shared_ptr<bam_hdr_t>& header) :
    m_header { header }
  {}

  /**
   * @brief creates a deep copy of a SamHeader object
   *
   * @note the copy will have exclusive ownership over the newly-allocated htslib memory
   */
  SamHeader::SamHeader(const SamHeader& other) :
    m_header { utils::make_shared_sam_header(utils::sam_header_deep_copy(other.m_header.get())) }
  {}

  /**
   * @brief moves a SamHeader object, transferring ownership of the underlying htslib memory
   */
  SamHeader::SamHeader(SamHeader&& other) noexcept :
    m_header { move(other.m_header) }
  {}

  /**
   * @brief creates a deep copy of a SamHeader object
   *
   * @note the copy will have exclusive ownership over the newly-allocated htslib memory
   */
  SamHeader& SamHeader::operator=(const SamHeader& other) {
    if ( &other == this )  
      return *this;
    m_header = utils::make_shared_sam_header(utils::sam_header_deep_copy(other.m_header.get())); ///< shared_ptr assignment will take care of deallocating old sam record if necessary
    return *this;
  }

  SamHeader& SamHeader::operator=(SamHeader&& other) noexcept {
    if (&other != this)
      m_header = move(other.m_header);
    return *this;
  }

  /**
   * @brief Returns the length of the given sequence as stored in the @SQ tag in the BAM header, or 0 if the sequence
   * name is not found.
   */
  uint32_t SamHeader::sequence_length(const std::string& sequence_name) const {
	  auto c = sequence_name.c_str();
	  for (int i = 0; i < m_header->n_targets; i++) {
		  if (strcmp(c,m_header->target_name[i]) == 0) {
			  return m_header->target_len[i];
		  }
	  }
	  return 0;
  }

}
