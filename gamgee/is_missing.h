#ifndef __gamgee__is_missing__
#define __gamgee__is_missing__

#include <cmath>

namespace gamgee {

bool is_missing(const float value) { return isnan(value) == 1; } ///< @brief checks whether or not a float value returned by a Variant field is missing

}

#endif
