#include "sam_utils.h"
#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>

namespace gamgee {

namespace utils {

/*
bool sam_index_exists(const std::string& filename) {
  // copy/pasted from http://www.cplusplus.com/forum/general/1796/
  const std::ifstream ifile(filename.c_str());
  return (bool)(ifile); // reinterpret_cast was not allowed by compiler, and link above said this should work... so... feel free to reimplement this function.
}

std::string sam_index_find(const std::string& filename) {
  // Check for extensions copied from sam_index_load in sam.h
  const std::vector<std::string> sam_index_extensions = {".csi", ".bai", ".crai"};
  for (const auto extension: sam_index_extensions) {
    const auto extension_appended = filename + extension;
    if (sam_index_exists(extension_appended)) {
      return extension_appended;
    }
  }

  // Check for the compact version of just swapping the last letter of the file with i.
  auto i_swapped = filename;
  i_swapped.replace(i_swapped.size()-1, 1, "i");
  if (sam_index_exists(i_swapped)) {
    return i_swapped;
  }

  throw std::invalid_argument("unable to locate index for: " + filename);
}
*/

}

}
