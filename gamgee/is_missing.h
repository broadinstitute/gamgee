#ifndef gamgee__is_missing__guard
#define gamgee__is_missing__guard

#include "sam_tag.h"

#include <cmath>

namespace gamgee {

bool is_missing(const float value);

template<class T>
bool is_missing(const SamTag<T>& tag) { return !tag.is_present(); } ///< @brief checks whether or not a SamTag is missing

}

#endif // gamgee__is_missing__guard
