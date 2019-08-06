#ifndef gamgee__header_line__guard
#define gamgee__header_line__guard

#include<string>
#include<unordered_map>

namespace gamgee {


/**
 * @brief Helper struct to hold header line from a sam header.
 */
struct SamHeaderLine {
  static constexpr char HD_LINE_CODE [] = "@HD";

  enum class SortingOrder {NOT_SET, UNKNOWN, UNSORTED, QUERYNAME, COORDINAMTE};
  enum class Grouping {NOT_SET, NONE, QUERY, REFERENCE};

  std::string version;
  SortingOrder sorting_order;
  Grouping grouping;
  std::unordered_map<std::string, std::string> extra_fields;

  SamHeaderLine() :
    version { DEFAULT_VERSION },
    sorting_order {SortingOrder::NOT_SET},
    grouping {Grouping::NOT_SET} {};
  SamHeaderLine(const SamHeaderLine& other) = default;
  SamHeaderLine(const std::string& header_line);

  bool operator==(const SamHeaderLine& header_line) const {
    bool result = ( (version == header_line.version) &&
                    (sorting_order == header_line.sorting_order) &&
                    (grouping == header_line.grouping));

    for (auto it = extra_fields.cbegin(); it != extra_fields.cend(); ++it) {
      result = result &&
               ( header_line.extra_fields.find(it->first) != header_line.extra_fields.cend() ) &&
               ( it->second == header_line.extra_fields.at(it->first) );
      if ( !result ) return result;
    }

    return result;
  }
  friend std::ostream& operator<<(std::ostream& out,
                                  const SamHeaderLine& header_line);

 private:
  static constexpr char DEFAULT_VERSION [] = "1.4";

};


}



#endif //gamgee__comment_guard
