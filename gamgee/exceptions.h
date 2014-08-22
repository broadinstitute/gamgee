#ifndef gamgee__exceptions__guard
#define gamgee__exceptions__guard

#include <boost/format.hpp>

#include <exception>
#include<string>

namespace gamgee {

/**
 * @brief an exception class for the case where a single input is required, but more is provided
 */
class SingleInputException : public std::runtime_error {
 public:
  SingleInputException(const std::string& vector_name, const size_t size) :
    std::runtime_error{(boost::format("Error: single input required, but vector %s has size %d") % vector_name % size).str()} { }
};

} // end of namespace gamgee

#endif // end of gamgee__exceptions__guard

