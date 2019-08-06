#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <map>

#include "header_line.h"
#include "../utils/utils.h"

using namespace std;
using namespace boost;

namespace gamgee {
  constexpr char SamHeaderLine::HD_LINE_CODE [];
  constexpr char SamHeaderLine::DEFAULT_VERSION [];

  constexpr char VERSION_TAG [] = "VN";
  constexpr char SORTING_TAG [] = "SO";
  constexpr char GROUPING_TAG [] = "GO";

  /**
   * @brief Creates and populates a SamHeaderLine object from a sam header line.
   *
   * @note The top header line in a SAM header starts with "@HD".
   * Example: @HD    SO:coordinate.
   */
  SamHeaderLine::SamHeaderLine(const std::string& header_line) {
    static constexpr auto CHARACTERS_PER_TAG = 2;
    static constexpr auto HD_TAG_LENGTH = 3;
    static const auto TAB_SEPARATOR = char_separator<char>("\t");

    // Init SO and GO fields.
    sorting_order = SortingOrder::NOT_SET;
    grouping = Grouping::NOT_SET;
    auto fields = header_line.substr(HD_TAG_LENGTH+1);  //skip the @HD tag
    auto tokens = tokenizer<char_separator<char>>(fields, TAB_SEPARATOR);
    for ( const auto& token : tokens) {
      auto tag = token.substr(0, CHARACTERS_PER_TAG);
      //we skip 3 characters -- 2 for the tag, 1 for the colon
      auto value = token.substr(CHARACTERS_PER_TAG+1, token.length() - CHARACTERS_PER_TAG - 1);
      if ( utils::starts_with(token, VERSION_TAG) ) version = value;
      else if ( utils::starts_with(token, SORTING_TAG) ) {
        if      (iequals(value, "unknown"))    sorting_order = SortingOrder::UNKNOWN;
        else if (iequals(value, "unsorted"))   sorting_order = SortingOrder::UNSORTED;
        else if (iequals(value, "queryname"))  sorting_order = SortingOrder::QUERYNAME;
        else if (iequals(value, "coordinate")) sorting_order = SortingOrder::COORDINAMTE;
        else throw logic_error("Unsupported sorting order (SO) in header line.");
      }
      else if ( utils::starts_with(token, GROUPING_TAG) ) {
        if      (iequals(value, "none"))      grouping = Grouping::NONE;
        else if (iequals(value, "query"))     grouping = Grouping::QUERY;
        else if (iequals(value, "reference")) grouping = Grouping::REFERENCE;
        else throw logic_error("Unsupported grouping (GO) in header line.");
      }
      else {
        extra_fields[tag] = value;
      }
    }
  }

  /**
   * @brief Overloads the << operator for SamHeaderLine object.
   *
   * @note A header line in a SAM header starts with "@HD". This function does not
   * check if header line string is valid.
   * Ex: @CO    My comment.
   */
  std::ostream& operator<<(std::ostream& out,
                           const SamHeaderLine& header_line) {

    out << SamHeaderLine::HD_LINE_CODE;

    if (!header_line.version.empty()) out << "\t" << VERSION_TAG << ":" << header_line.version;
    if (header_line.sorting_order != SamHeaderLine::SortingOrder::NOT_SET) {
      out << "\t" << SORTING_TAG << ":";
      switch(header_line.sorting_order) {
        case SamHeaderLine::SortingOrder::UNKNOWN:     out << "unknown"; break;
        case SamHeaderLine::SortingOrder::UNSORTED:    out << "unsorted"; break;
        case SamHeaderLine::SortingOrder::QUERYNAME:   out << "queryname"; break;
        case SamHeaderLine::SortingOrder::COORDINAMTE: out << "coordinate"; break;
        case SamHeaderLine::SortingOrder::NOT_SET: break;
      }
    }
    if (header_line.grouping != SamHeaderLine::Grouping::NOT_SET) {
      out << "\t" << GROUPING_TAG << ":";
      switch(header_line.grouping) {
        case SamHeaderLine::Grouping::NONE:      out << "none"; break;
        case SamHeaderLine::Grouping::QUERY:     out << "query"; break;
        case SamHeaderLine::Grouping::REFERENCE: out << "reference"; break;
        case SamHeaderLine::Grouping::NOT_SET: break;
      }
    }
    // For consistency, add extra fields in alphabetical order.
    std::map<std::string, std::string> ordered_fields(
        header_line.extra_fields.cbegin(), header_line.extra_fields.cend());
    for (const auto& field : ordered_fields) {
      out << "\t"  << field.first  << ":"  << field.second;
    }
    return out;
  }

}
