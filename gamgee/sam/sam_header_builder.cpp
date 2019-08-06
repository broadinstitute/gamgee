
#include "sam_header_builder.h"

#include "../utils/hts_memory.h"

#include <algorithm>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>

using namespace std;

namespace gamgee {

  /**
   * @brief create a SamHeader from scratch.
   *
   * @note users of this constructor will need to fill in sufficient parts of the header
   *       to produce a valid header.
   */
  SamHeaderBuilder::SamHeaderBuilder(const bool validate_on_build) :
    sequences_info {},
    read_groups {},
    programs {},
    comments {},
    m_validate_on_build { validate_on_build }
  {}

  /**
   * @brief create a SamHeader from an existing header.
   *
   * @note the data in the existing header is copied into the builder so that future changes
   *       to the header will not affect the builder.
   */
  SamHeaderBuilder::SamHeaderBuilder(const SamHeader& sam_header, bool validate_on_build) :
    header_line { sam_header.header_line() },
    read_groups { sam_header.read_groups() },
    programs { sam_header.programs() },
    comments { sam_header.comments()},
    m_validate_on_build { validate_on_build }
  {
    sequences_info.resize(sam_header.n_sequences());
    for (int i = 0; i != sam_header.n_sequences(); ++i) {
      std::string seq_name = sam_header.sequence_name(i);
      uint32_t seq_len = sam_header.sequence_length(i);
      sequences_info[i] = std::make_pair(seq_name, seq_len);
    }
  }

  /**
   * @brief Set header line information.
   *
   * @note: set_header_line() makes a copy of input data, so its future changes will not
   * affect the builder.
   */
  SamHeaderBuilder& SamHeaderBuilder::set_header_line(
      const SamHeaderLine& header_line_) {
    header_line = SamHeaderLine(header_line_);
    return *this;
  }

  /**
   * @brief Set reference sequences information.
   *
   * @note: set_seqs_info() makes a copy of input data, so its future changes will not
   * affect the builder.
   */
  SamHeaderBuilder& SamHeaderBuilder::set_seqs_info(
      const std::vector<std::pair<std::string, uint32_t>>& seqs_info) {
    sequences_info = std::vector<std::pair<std::string, uint32_t>>();
    for (const auto& seq_info : seqs_info) append_seq_info(seq_info.first, seq_info.second);
    return *this;
  }

  /**
   * @brief Append reference sequence information.
   *
   * @note: append_seq_info() makes a copy of input data, so its future changes will not
   * affect the builder.
   */
  SamHeaderBuilder& SamHeaderBuilder::append_seq_info(string seq_name, uint32_t seq_len) {
    sequences_info.push_back(std::make_pair(seq_name, seq_len));
    return *this;
  }

  /**
   * @brief Set programs info.
   *
   * @note: seq_programs() makes a copy of input data, so its future changes will not
   * affect the builder.
   */
  SamHeaderBuilder& SamHeaderBuilder::set_programs(const std::vector<Program>& programs_) {
    programs = std::vector<Program>();
    for (const auto& program : programs_) append_program(program);
    return *this;
  }

  /**
   * @brief Append a single program info.
   *
   * @note: append_program() makes a copy of input data, so its future changes will not
   * affect the builder.
   */
  SamHeaderBuilder& SamHeaderBuilder::append_program(const Program& program) {
    programs.push_back(Program(program));
    return *this;
  }

  /**
   * @brief Set read groups info.
   *
   * @note: set_read_groups() makes a copy of input data, so its future changes will not
   * affect the builder.
   */
  SamHeaderBuilder& SamHeaderBuilder::set_read_groups(const std::vector<ReadGroup>& read_groups_) {
    read_groups = std::vector<ReadGroup>();
    for (const auto& read_group : read_groups_) append_read_group(read_group);
    return *this;
  }

  /**
   * @brief Append a single read group info.
   *
   * @note: append_read_group() makes a copy of input data, so its future changes will not
   * affect the builder.
   */
  SamHeaderBuilder& SamHeaderBuilder::append_read_group(const ReadGroup& read_group) {
    read_groups.push_back(ReadGroup(read_group));
    return *this;
  }

  /**
   * @brief Set comments.
   *
   * @note: set_comments() makes a copy of input data, so its future changes will not
   * affect the builder.
   */
  SamHeaderBuilder& SamHeaderBuilder::set_header_comments(const std::vector<SamHeaderComment>& comments_) {
    comments = std::vector<SamHeaderComment>();
    for (const auto& comment : comments_) append_header_comment(comment);
    return *this;
  }

  /**
   * @brief Append a single comment.
   *
   * @note: append_comment() makes a copy of input data, so its future changes will not
   * affect the builder.
   */
  SamHeaderBuilder& SamHeaderBuilder::append_header_comment(const SamHeaderComment& comment) {
    comments.push_back(SamHeaderComment(comment));
    return *this;
  }

  /**
   * @brief build a SamHeader record using the current state of the builder
   *
   * @note this version of build() can be called repeatedly to build multiple SamHeader objects.
   */
  SamHeader SamHeaderBuilder::build() const {
    // Perform validation before building, if requested
    if ( m_validate_on_build )
      validate();

    std::stringstream text_stream;
    text_stream << header_line << "\n";
    std::copy(read_groups.cbegin(), read_groups.cend(), std::ostream_iterator<ReadGroup>(text_stream, "\n"));
    std::copy(programs.cbegin(), programs.cend(), std::ostream_iterator<Program>(text_stream, "\n"));
    std::copy(comments.cbegin(), comments.cend(), std::ostream_iterator<SamHeaderComment>(text_stream, "\n"));

    bam_hdr_t* bam_header = bam_hdr_init();

    bam_header->n_targets = sequences_info.size();
    bam_header->target_name = (char**) malloc(sequences_info.size() * sizeof(char*));
    bam_header->target_len = (uint32_t*) malloc(sequences_info.size() * sizeof(uint32_t));

    int contig_i = 0;
    for (const auto& sequence_info : sequences_info) {
      bam_header->target_name[contig_i] = (char*) malloc((sequence_info.first.size() + 1) * sizeof(char));
      memcpy(bam_header->target_name[contig_i], sequence_info.first.c_str(),
             sequence_info.first.size());
      bam_header->target_name[contig_i][sequence_info.first.size()] = '\0';
      bam_header->target_len[contig_i] = sequence_info.second;
      contig_i++;
    }

    auto text = text_stream.str();
    bam_header->l_text = text.length();
    bam_header->text = (char*) malloc((text.length() + 1) * sizeof(char));
    memcpy(bam_header->text, text.c_str(), text.length());
    bam_header->text[text.length()] = '\0';

    gamgee::SamHeader sam_header(
        gamgee::utils::make_shared_sam_header(bam_header));

    return sam_header;
  }

  /**
   * @brief performs pre-build validation of the state of the SamHeader record under construction
   */
  void SamHeaderBuilder::validate() const {
    // For any reference squence info, sequence name should not be empty.
    for (const auto& seq : sequences_info) {
      if (seq.first.empty())
        throw logic_error("Sam header reference name is an empty string.");
    }

    // For any read group info, ID is a required field.
    for (const auto& read_group: read_groups) {
      if (read_group.id.empty())
        throw logic_error("Missing required ID field for sam header read group info.");
    }

    // For any read program , ID is a required field.
    for (const auto& program: programs) {
      if (program.id.empty())
        throw logic_error("Missing required ID field for sam header program info.");
    }

    // Each comment string should be a single-line.
    for (const auto& comment: comments) {
      if (comment.comment.find('\n') != std::string::npos)
        throw logic_error("Sam header comment string should not be multi-line.");
    }

    // TODO: add additional validation for fileds tag/value regex match.
  }

}
