#include <boost/tokenizer.hpp>
#include <map>

#include "program.h"
#include "../utils/utils.h"

using namespace std;
using namespace boost;

namespace gamgee {
  constexpr char Program::PG_LINE_CODE [];

  constexpr char ID_TAG [] = "ID";
  constexpr char NAME_TAG [] = "PN";
  constexpr char COMMAND_LINE_TAG [] = "CL";
  constexpr char VERSION_TAG [] = "VN";

  /**
   * @brief creates and populates a program object from a sam header line
   *
   * @note A program line in a SAM header starts with "@PG" and then has tag/value
   * pairs separated by tab.  The tags are all length-2 strings and the tags are
   * separated from values by colons.
   * Ex: @PG     ID:bwa  PN:bwa  VN:0.7.12-r1039 CL:bwa mem -t 12 index.fa read1.fq read2.fq
   */
  Program::Program(const std::string& header_line) {
    static constexpr auto CHARACTERS_PER_TAG = 2;
    static constexpr auto PG_TAG_LENGTH = 3;
    static const  auto TAB_SEPARATOR = char_separator<char>("\t");

    auto fields = header_line.substr(PG_TAG_LENGTH+1);  //skip the @PG tag
    auto tokens = tokenizer<char_separator<char>>(fields, TAB_SEPARATOR);
    for ( const auto& token : tokens) {
      auto tag = token.substr(0, CHARACTERS_PER_TAG);
      //we skip 3 characters -- 2 for the tag, 1 for the colon
      auto value = token.substr(CHARACTERS_PER_TAG+1, token.length() - CHARACTERS_PER_TAG - 1);
      if      ( utils::starts_with(token, ID_TAG) )            id = value;
      else if ( utils::starts_with(token, NAME_TAG) )          name = value;
      else if ( utils::starts_with(token, VERSION_TAG) )       version = value;
      else if ( utils::starts_with(token, COMMAND_LINE_TAG) )  command_line = value;
      else extra_fields[tag] = value;
    }
  }

  /**
   * @brief Overloads the << operator for Program object.
   *
   * @note A program info line in a SAM header starts with "@PG" and then has tag/value
   * pairs separated by tab. The tags are all length-2 strings and the tags are
   * separated from values by colons. This function does not check if program line string is valid.
   * Ex: @PG     ID:bwa  PN:bwa  VN:0.7.12-r1039 CL:bwa mem -t 12 index.fa read1.fq read2.fq
   */
  std::ostream& operator<<(std::ostream& out,
                           const Program& program) {
    out << Program::PG_LINE_CODE;

    if (!program.id.empty())             out << "\t" << ID_TAG << ":" << program.id;
    if (!program.name.empty())           out << "\t" << NAME_TAG << ":" << program.name;
    if (!program.command_line.empty())   out << "\t" << COMMAND_LINE_TAG << ":" << program.command_line;
    if (!program.version.empty())        out << "\t" << VERSION_TAG << ":" << program.version;
    // For consistency, add extra fields in alphabetical order.
    std::map<std::string, std::string> ordered_fields(
        program.extra_fields.cbegin(), program.extra_fields.cend());
    for (const auto& field : ordered_fields) {
      out << "\t" << field.first << ":" << field.second;
    }
    return out;
  }

}
