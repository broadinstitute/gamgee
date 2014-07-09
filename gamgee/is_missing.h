#ifndef __gamgee__is_missing__
#define __gamgee__is_missing__

#include "sam_tag.h"

#include <cmath>

namespace gamgee {

bool is_missing(const float value);

template<class T>
bool is_missing(const SamTag<T>& tag) { return !tag.is_present(); } ///< @brief checks whether or not a SamTag is missing

}

#endif
