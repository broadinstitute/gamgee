#ifndef gamgee__zip__guard
#define gamgee__zip__guard

#include <boost/format.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

namespace gamgee {
namespace utils {

/**
 * @brief utility method to zip iterators together with simpler syntax than boost
 *
 * This is a wrapper over boost's zip_iterator interface to simplify the usage of zip
 * iterators especially in for each loops. This function enables the following syntax:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (const auto tup : zip(a, b, c, d) {
 *   ... // use tup values as a boost::tuple
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * for more details look at boost's zip_iterator documentation.
 */
template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(containers.begin()...))>>
{
    auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
    auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
    return boost::make_iterator_range(zip_begin, zip_end);
}

} // namespace utils
} // namespace gamgee

#endif /* gamgee__zip__guard */
