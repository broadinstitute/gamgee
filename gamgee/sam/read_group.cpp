#include <sstream>

#include <boost/tokenizer.hpp>

#include "read_group.h"

using namespace std;
using namespace boost;

namespace gamgee {


  /**
   * @brief creates and populates a read group object from a sam header line
   *
   * @note A read group line in a SAM header starts with "@RG" and then has tag/value
   * pairs separated by whitespace.  The tags are all length-2 strings and the tags are
   * separated from values by colons.  Ex: @RG   ID:12345.67   SM:11111  PL:Illumina
   */
  ReadGroup::ReadGroup(const std::string& header_line) {
    const static auto CHARACTERS_PER_TAG = 2;
    const static auto RG_TAG_LENGTH = 3;
    const static auto TAB_SEPARATOR = char_separator<char>("\t");

    auto fields = header_line.substr(RG_TAG_LENGTH+1);  //skip the @RG tag
    auto tokens = tokenizer<char_separator<char>>(fields, TAB_SEPARATOR);
    for ( const auto& token : tokens) {
      //we skip 3 characters -- 2 for the tag, 1 for the colon
      auto value = token.substr(CHARACTERS_PER_TAG+1, token.length() - CHARACTERS_PER_TAG - 1);
      if      ( starts_with(token, ID_TAG) )                  id = value;
      else if ( starts_with(token, SAMPLE_TAG) )              sample = value;
      else if ( starts_with(token, CENTER_TAG) )              center = value;
      else if ( starts_with(token, DESCRIPTION_TAG) )         description = value;
      else if ( starts_with(token, DATE_TIME_TAG) )           date_time = value;
      else if ( starts_with(token, CENTER_TAG) )              center = value;
      else if ( starts_with(token, FLOW_ORDER_TAG) )          flow_order = value;
      else if ( starts_with(token, KEY_SEQUENCE_TAG) )        key_sequence = value;
      else if ( starts_with(token, LIBRARY_TAG) )             library = value;
      else if ( starts_with(token, PROGRAMS_TAG) )            programs = value;
      else if ( starts_with(token, MEDIAN_INSERT_SIZE_TAG) )  median_insert_size = value;
      else if ( starts_with(token, PLATFORM_TAG) )            platform = value;
      else if ( starts_with(token, PLATFORM_UNIT_TAG) )       platform_unit = value;
    }
  }

  const char ReadGroup::ID_TAG [3] = "ID";
  const char ReadGroup::CENTER_TAG [3] = "CN";
  const char ReadGroup::DESCRIPTION_TAG [3] = "DS";
  const char ReadGroup::DATE_TIME_TAG [3] = "DT";
  const char ReadGroup::FLOW_ORDER_TAG [3] = "FO";
  const char ReadGroup::KEY_SEQUENCE_TAG [3] = "KS";
  const char ReadGroup::LIBRARY_TAG [3] = "LB";
  const char ReadGroup::PROGRAMS_TAG [3] = "PG";
  const char ReadGroup::MEDIAN_INSERT_SIZE_TAG [3] = "PI";
  const char ReadGroup::PLATFORM_TAG [3] = "PL";
  const char ReadGroup::PLATFORM_UNIT_TAG [3] = "PU";
  const char ReadGroup::SAMPLE_TAG [3] = "SM";
}
