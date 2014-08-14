#ifndef gamgee__index_file_utils__guard
#define gamgee__index_file_utils__guard

#include <string>
#include <vector>

namespace gamgee {

namespace utils {

namespace contig_values {
  const auto entire_file = std::string{"."}; // Copied from strcmp line in hts_itr_querys()
  const auto entire_file_list = std::vector<std::string>{entire_file};
  const auto unmapped_contig = std::string{"*"}; // Copied from strcmp line in hts_itr_querys(). Very, very different from unmapped flagged reads.
}

/*
bool sam_index_exists(const std::string& filename);
std::string sam_index_find(const std::string& filename);
*/

}

}

#endif
