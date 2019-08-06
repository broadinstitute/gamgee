#ifndef gamgee__sam_header_builder__guard
#define gamgee__sam_header_builder__guard

#include "header_line.h"
#include "read_group.h"
#include "program.h"
#include "header_comment.h"
#include "sam_header.h"

#include <string>
#include <utility>
#include <vector>

namespace gamgee {

/**
 * @brief class to build SamHeader objects from existing data or from scratch.
 *
 * To use, either start from scratch:
 *
 * auto builder = SamHeaderBuilder{};
 *
 * or start with an exisiting header object:
 *
 * auto builder = SamHeaderBuilder{existing_header};
 *
 * Then make the desired modifications and build repeatedly.
 *
 * You must provide a value for all essential data fields before building if you
 * build from scratch, unless you disable validation (not recommended).
 *
 * The builder always copies the data it's given in the set_* methods, for safety reasons.
 * Altering a header after passing some of its data into a builder will not affect the
 * state of the builder.
 *
 * Before building, the builder performs a validation step to ensure that the read it will
 * create is logically consistent. This validation step can be disabled during construction
 * of the builder, but disabling validation is dangerous and not recommended.
 */
class SamHeaderBuilder {
 public:
  explicit SamHeaderBuilder(const bool validate_on_build = true);                           ///< @brief create a SamHeader from scratch.
  explicit SamHeaderBuilder(const SamHeader& header, const bool validate_on_build = true);  ///< @brief create a SamHeader with an existing header.

  // SamHeaderBuilders are moveable but not copyable, and use default destruction
  SamHeaderBuilder(SamHeaderBuilder&& other) = default;
  SamHeaderBuilder& operator=(SamHeaderBuilder&& other) = default;
  SamHeaderBuilder(const SamHeaderBuilder& other) = delete;
  SamHeaderBuilder& operator=(const SamHeaderBuilder& other) = delete;
  ~SamHeaderBuilder() = default;

  // Setters
  SamHeaderBuilder& set_header_line(const SamHeaderLine& header_line);
  SamHeaderBuilder& set_seqs_info(const std::vector<std::pair<std::string, uint32_t>>& seqs_info);
  SamHeaderBuilder& set_programs(const std::vector<Program>& programs);
  SamHeaderBuilder& set_read_groups(const std::vector<ReadGroup>& read_groups);
  SamHeaderBuilder& set_header_comments(const std::vector<SamHeaderComment>& comments);

  SamHeaderBuilder& append_seq_info(std::string seq_name, uint32_t seq_len);
  SamHeaderBuilder& append_program(const Program& program);
  SamHeaderBuilder& append_read_group(const ReadGroup& read_group);
  SamHeaderBuilder& append_header_comment(const SamHeaderComment& comment);

  SamHeader build() const;    ///< build a SamHeader (can be called repeatedly)

 private:
  std::vector<std::pair<std::string, uint32_t>> sequences_info;      ///< list of refernce sequences info (@SQ).
  SamHeaderLine header_line;                        ///< header line info (@HD) to be added to header text.
  std::vector<ReadGroup> read_groups;            ///< list of read groups (@RG) to be added to header text.
  std::vector<Program> programs;                 ///< list of programs (@PG) to be added to header text.
  std::vector<SamHeaderComment> comments;                 ///< list of comment lines (@CO) to be added to header text.
  bool m_validate_on_build;                      ///< should we validate the state of the SamHeader record at build time?

  void validate() const;                    ///< @brief performs pre-build validation of the state of the SamHeader record under construction.
};

}


#endif /* gamgee__sam_header_builder__guard */
