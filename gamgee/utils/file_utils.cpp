#include "file_utils.h"

#include <memory>
#include <fstream>
#include <string>

using namespace std;

namespace gamgee {
namespace utils {

std::shared_ptr<std::ifstream> make_shared_ifstream(std::ifstream* ifstream_ptr) {
  return shared_ptr<std::ifstream>(ifstream_ptr, IFStreamDeleter());
}

std::shared_ptr<std::ifstream> make_shared_ifstream(std::string filename) {
  return make_shared_ifstream(new std::ifstream{filename});
}

}
}
