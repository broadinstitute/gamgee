#ifndef gamgee__file_utils__guard
#define gamgee__file_utils__guard

#include <memory>
#include <fstream>
#include <string>

namespace gamgee {
namespace utils {

/**
 * @brief a functor object to delete an ifstream
 */
struct IFStreamDeleter {
  void operator()(std::ifstream* p) const { p->close(); }
};

/**
  * @brief wraps a pre-allocated ifstream in a shared_ptr with correct deleter
  * @param ifstream_ptr an ifstream raw file pointer
  */
std::shared_ptr<std::ifstream> make_shared_ifstream(std::ifstream* ifstream_ptr);

/**
  * @brief wraps an input file in a shared_ptr to an ifstream with correct deleter
  * @param filename the input filename
  */
std::shared_ptr<std::ifstream> make_shared_ifstream(std::string filename);

}
}

#endif /* gamgee__file_utils__guard */
