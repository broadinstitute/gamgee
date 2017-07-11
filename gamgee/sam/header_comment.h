#ifndef gamgee__header_comment__guard
#define gamgee__header_comment__guard

#include<string>

namespace gamgee {


/**
 * @brief Helper struct to hold one comment line from a sam header.
 */
struct SamHeaderComment {
  static constexpr char CO_LINE_CODE [] = "@CO";
  std::string comment;

  SamHeaderComment() = default;
  SamHeaderComment(const SamHeaderComment& other) = default;
  SamHeaderComment(const std::string& header_line);

  bool operator==(const SamHeaderComment& header_comment) const {
    return(comment == header_comment.comment);
  }
  friend std::ostream& operator<<(std::ostream& out,
                                  const SamHeaderComment& header_comment);
};


}



#endif //gamgee__header_comment_guard
