#ifndef gamgee__program__guard
#define gamgee__program__guard

#include<string>
#include<unordered_map>

namespace gamgee {


/**
 * @brief Helper struct to hold one program record from a sam file header
 *
 * A program, which is contained in a line starting with "@PG" in a sam file header,
 * contains information about the programs run on the dataset with the parameters used.
 */
struct Program {
  static constexpr const char PG_LINE_CODE [] = "@PG";

  std::string id;
  std::string name;
  std::string version;
  std::string command_line;
  std::unordered_map<std::string, std::string> extra_fields;

  Program() = default;
  Program(const Program& other) = default;
  Program(const std::string& program);

  bool operator==(const Program& program) const {
    bool result = ( (id == program.id) &&
                    (name == program.name) &&
                    (version == program.version) &&
                    (command_line == program.command_line));
    for (auto it = extra_fields.cbegin(); it != extra_fields.cend(); ++it) {
      result = result &&
               ( program.extra_fields.find(it->first) != program.extra_fields.cend() ) &&
               ( it->second == program.extra_fields.at(it->first) );
      if ( !result ) return result;
    }

    return result;
  }
  friend std::ostream& operator<<(std::ostream& out,
                                  const Program& program);

};


}



#endif //gamgee__program__guard
