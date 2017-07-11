
#include <sstream>

#include <boost/tokenizer.hpp>

#include "header_comment.h"

using namespace std;
using namespace boost;

namespace gamgee {
  constexpr char SamHeaderComment::CO_LINE_CODE [];

  /**
   * @brief creates and populates a SamHeaderComment object from a sam header line.
   *
   * @note A comment line in a SAM header starts with "@CO".
   * Example: @CO    My comment.
   */
  SamHeaderComment::SamHeaderComment(const std::string& header_line) {
    static constexpr auto CO_TAG_LENGTH = 3;
    static const auto TAB_SEPARATOR = char_separator<char>("\t");
    comment = header_line.substr(CO_TAG_LENGTH+1);  //skip the @CO\t tag
  }

  /**
   * @brief Overloads the << operator for SamHeaderComment object.
   *
   * @note A comment line in a SAM header starts with "@CO". This function does
   * not check if comment line string is valid.
   * Ex: @CO    My comment.
   */
  std::ostream& operator<<(std::ostream& out,
                           const SamHeaderComment& header_comment) {
    out << SamHeaderComment::CO_LINE_CODE << '\t' << header_comment.comment;
    return out;
  }

}
