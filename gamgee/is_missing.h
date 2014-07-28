#ifndef gamgee__is_missing__guard
#define gamgee__is_missing__guard

#include "sam_tag.h"

#include <string>
#include <cmath>

namespace gamgee {

inline bool is_missing(const float value) { return isnan(value) == 1; }                    ///< @brief checks whether or not a float value is representing missing data. This works for all floats in Gamgee.
inline bool is_missing(const std::string& value) { return value.empty() || value == "."; } ///< @brief checks whether or not a string value returned by a Variant field is empty or equal to ".". This works for all strings in Gamgee

template<class T>
bool is_missing(const SamTag<T>& tag) { return !tag.is_present(); }                        ///< @brief checks whether or not a SamTag is missing

}

#endif // gamgee__is_missing__guard
