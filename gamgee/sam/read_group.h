#ifndef gamgee__read_group__guard
#define gamgee__read_group__guard

#include<string>
#include<unordered_map>

namespace gamgee {


/**
 * @brief Helper struct to hold one read group record from a sam file header
 *
 * A read group, which is contained in a line starting with "@RG" in a sam file header,
 * contains information about the seuencing process: type of machine, time, center etc.
 * Most important for variant calling and quality control is the sample ID -- each sample
 * may be split into several read groups.
 */
struct ReadGroup {
  static constexpr const char RG_LINE_CODE [] = "@RG";

  std::string id;
  std::string center;
  std::string description;
  std::string date_time;
  std::string flow_order;
  std::string key_sequence;
  std::string library;
  std::string programs;
  std::string median_insert_size;
  std::string platform;
  std::string platform_model;
  std::string platform_unit;
  std::string sample;
  std::unordered_map<std::string, std::string> extra_fields;

  ReadGroup() = default;
  ReadGroup(const ReadGroup& other) = default;
  ReadGroup(const std::string& header_line);

  bool operator==(const ReadGroup& read_group) const {
    bool result = ( (id == read_group.id) &&
                    (center == read_group.center) &&
                    (description == read_group.description) &&
                    (date_time == read_group.date_time) &&
                    (flow_order == read_group.flow_order) &&
                    (key_sequence == read_group.key_sequence) &&
                    (library == read_group.library) &&
                    (programs == read_group.programs) &&
                    (median_insert_size == read_group.median_insert_size) &&
                    (platform == read_group.platform) &&
                    (platform_model == read_group.platform_model) &&
                    (platform_unit == read_group.platform_unit) &&
                    (sample == read_group.sample));
    for (auto it = extra_fields.cbegin(); it != extra_fields.cend(); ++it) {
      result = result &&
               ( read_group.extra_fields.find(it->first) != read_group.extra_fields.cend() ) &&
               ( it->second == read_group.extra_fields.at(it->first));
      if ( !result ) return result;
    }

    return(result);
  }
  friend std::ostream& operator<<(std::ostream& out,
                                  const ReadGroup& read_group);

};


}



#endif //gamgee__read_group__guard
