#include <sstream>

#include <boost/tokenizer.hpp>
#include <map>

#include "read_group.h"
#include "../utils/utils.h"

using namespace std;
using namespace boost;

namespace gamgee {
  constexpr char ReadGroup::RG_LINE_CODE [];

  constexpr char ID_TAG [] = "ID";
  constexpr char CENTER_TAG [] = "CN";
  constexpr char DESCRIPTION_TAG [] = "DS";
  constexpr char DATE_TIME_TAG [] = "DT";
  constexpr char FLOW_ORDER_TAG [] = "FO";
  constexpr char KEY_SEQUENCE_TAG [] = "KS";
  constexpr char LIBRARY_TAG [] = "LB";
  constexpr char PROGRAMS_TAG [] = "PG";
  constexpr char MEDIAN_INSERT_SIZE_TAG [] = "PI";
  constexpr char PLATFORM_TAG [] = "PL";
  constexpr char PLATFORM_MODEL_TAG [] = "PM";
  constexpr char PLATFORM_UNIT_TAG [] = "PU";
  constexpr char SAMPLE_TAG [] = "SM";

  /**
   * @brief creates and populates a read group object from a sam header line
   *
   * @note A read group line in a SAM header starts with "@RG" and then has tag/value
   * pairs separated by tab.  The tags are all length-2 strings and the tags are
   * separated from values by colons.  Ex: @RG   ID:12345.67   SM:11111  PL:Illumina
   */
  ReadGroup::ReadGroup(const std::string& header_line) {
    static constexpr auto CHARACTERS_PER_TAG = 2;
    static constexpr auto RG_TAG_LENGTH = 3;
    static const auto TAB_SEPARATOR = char_separator<char>("\t");

    auto fields = header_line.substr(RG_TAG_LENGTH+1);  //skip the @RG tag
    auto tokens = tokenizer<char_separator<char>>(fields, TAB_SEPARATOR);
    for ( const auto& token : tokens) {
      auto tag = token.substr(0, CHARACTERS_PER_TAG);
      //we skip 3 characters -- 2 for the tag, 1 for the colon
      auto value = token.substr(CHARACTERS_PER_TAG+1, token.length() - CHARACTERS_PER_TAG - 1);
      if      ( utils::starts_with(token, ID_TAG) )                  id = value;
      else if ( utils::starts_with(token, CENTER_TAG) )              center = value;
      else if ( utils::starts_with(token, DESCRIPTION_TAG) )         description = value;
      else if ( utils::starts_with(token, DATE_TIME_TAG) )           date_time = value;
      else if ( utils::starts_with(token, FLOW_ORDER_TAG) )          flow_order = value;
      else if ( utils::starts_with(token, KEY_SEQUENCE_TAG) )        key_sequence = value;
      else if ( utils::starts_with(token, LIBRARY_TAG) )             library = value;
      else if ( utils::starts_with(token, PROGRAMS_TAG) )            programs = value;
      else if ( utils::starts_with(token, MEDIAN_INSERT_SIZE_TAG) )  median_insert_size = value;
      else if ( utils::starts_with(token, PLATFORM_TAG) )            platform = value;
      else if ( utils::starts_with(token, PLATFORM_MODEL_TAG) )      platform_model = value;
      else if ( utils::starts_with(token, PLATFORM_UNIT_TAG) )       platform_unit = value;
      else if ( utils::starts_with(token, SAMPLE_TAG) )              sample = value;
      else extra_fields[tag] = value;
    }
  }

  /**
   * @brief Overloads the << operator for ReadGroup object.
   *
   * @note A read group line in a SAM header starts with "@RG" and then has tag/value
   * pairs separated by tab.  The tags are all length-2 strings and the tags are
   * separated from values by colons. This function does not check if read group line string
   * is valid.
   * Ex: @RG   ID:12345.67   SM:11111  PL:Illumina
   */
  std::ostream& operator<<(std::ostream& out,
                           const ReadGroup& read_group) {
    out << ReadGroup::RG_LINE_CODE;

    if (!read_group.id.empty())            out << "\t" << ID_TAG << ":" << read_group.id;
    if (!read_group.center.empty())        out << "\t" << CENTER_TAG << ":" << read_group.center;
    if (!read_group.description.empty())   out << "\t" << DESCRIPTION_TAG << ":" << read_group.description;
    if (!read_group.date_time.empty())     out << "\t" << DATE_TIME_TAG << ":" << read_group.date_time;
    if (!read_group.flow_order.empty())    out << "\t" << FLOW_ORDER_TAG << ":" << read_group.flow_order;
    if (!read_group.key_sequence.empty())  out << "\t" << KEY_SEQUENCE_TAG << ":" << read_group.key_sequence;
    if (!read_group.library.empty())       out << "\t" << LIBRARY_TAG << ":" << read_group.library;
    if (!read_group.programs.empty())      out << "\t" << PROGRAMS_TAG << ":" << read_group.programs;
    if (!read_group.median_insert_size.empty())  out << "\t" << MEDIAN_INSERT_SIZE_TAG << ":" << read_group.median_insert_size;
    if (!read_group.platform.empty())      out << "\t" << PLATFORM_TAG << ":" << read_group.platform;
    if (!read_group.platform_model.empty())      out << "\t" << PLATFORM_MODEL_TAG << ":" << read_group.platform_model;
    if (!read_group.platform_unit.empty()) out << "\t" << PLATFORM_UNIT_TAG << ":" << read_group.platform_unit;
    if (!read_group.sample.empty())        out << "\t" << SAMPLE_TAG << ":" << read_group.sample;
    // For consistency, add extra fields in alphabetical order.
    std::map<std::string, std::string> ordered_fields(
        read_group.extra_fields.cbegin(), read_group.extra_fields.cend());
    for (const auto& field : ordered_fields) {
      out << "\t" << field.first << ":" << field.second;
    }
    return out;
  }
}
