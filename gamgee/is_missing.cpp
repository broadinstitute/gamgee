#include "is_missing.h"

#include <cmath>

namespace gamgee {

bool is_missing(const float value) { return isnan(value) == 1; } ///< @brief checks whether or not a float value returned by a Variant field is missing

}
