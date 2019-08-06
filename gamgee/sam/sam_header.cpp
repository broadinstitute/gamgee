#include <boost/tokenizer.hpp>
#include "sam_header.h"

#include "../utils/hts_memory.h"
#include "../utils/utils.h"

using namespace std;
using namespace boost;

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

  /**
   * @brief Returns the length of the given sequence as stored in the \@SQ tag in the BAM header, or 0 if the sequence
   * name is not found.
   */
  uint32_t SamHeader::sequence_length(const std::string& sequence_name) const {
    auto c = sequence_name.c_str();
    for (int i = 0; i != m_header->n_targets; ++i) {
      if (strcmp(c,m_header->target_name[i]) == 0) {
        return m_header->target_len[i];
      }
    }
    return 0;
  }

  /**
   * @brief extracts SamHeaderLine object from a SAM header
   */
  SamHeaderLine SamHeader::header_line() const {
    auto result = SamHeaderLine();
    const auto text = header_text();

    // header line (@HD) should be the first line, if present.
    if (utils::starts_with(text, SamHeaderLine::HD_LINE_CODE) ) {
      auto hd_end = text.find('\n');
      auto hd_record = text.substr(0, hd_end);
      result = SamHeaderLine(hd_record);
    }
    return result;
  }

  /**
   * @brief extracts read group objects from a SAM header
   */
  vector<ReadGroup> SamHeader::read_groups() const {
    const static auto LINE_SEPARATOR = char_separator<char>("\n");
    auto result = vector<ReadGroup>();
    const auto text = header_text();

    const auto lines = tokenizer<char_separator<char>>(text, LINE_SEPARATOR);
    for (const auto& line : lines) {
      if ( utils::starts_with(line, ReadGroup::RG_LINE_CODE) ) {
        result.push_back(ReadGroup{line});
      }
    }
    return result;
  }

  /**
   * @brief extracts program objects from a SAM header
   */
  vector<Program> SamHeader::programs() const {
    const static auto LINE_SEPARATOR = char_separator<char>("\n");
    auto result = vector<Program>();
    const auto text = header_text();

    const auto lines = tokenizer<char_separator<char>>(text, LINE_SEPARATOR);
    for (auto& line : lines) {
      if ( utils::starts_with(line, Program::PG_LINE_CODE) ) {
        result.push_back(Program{line});
      }
    }
    return result;
  }

  /**
   * @brief extracts comments objects from a SAM header
   */
  vector<SamHeaderComment> SamHeader::comments() const {
    const static auto LINE_SEPARATOR = char_separator<char>("\n");
    auto result = vector<SamHeaderComment>();
    const auto text = header_text();

    const auto lines = tokenizer<char_separator<char>>(text, LINE_SEPARATOR);
    for (const auto& line : lines) {
      if ( utils::starts_with(line, SamHeaderComment::CO_LINE_CODE) ) {
        result.push_back(SamHeaderComment{line});
      }
    }
    return result;
  }
}
