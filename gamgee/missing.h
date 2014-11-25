#ifndef gamgee__missing__guard
#define gamgee__missing__guard

#include "sam/sam_tag.h"
#include "utils/utils.h"
#include "htslib/vcf.h"

#include <vector>
#include <string>
#include <cmath>

namespace gamgee {

namespace missing_values {
  constexpr auto int8 = bcf_int8_missing;                                                                                       ///< missing value for an int8
  constexpr auto int16 = bcf_int16_missing;                                                                                     ///< missing value for an int16
  constexpr auto int32 = bcf_int32_missing;                                                                                     ///< missing value for an int32
  constexpr auto string_empty = "";                                                                                             ///< empty string is a missing string
  constexpr auto string_dot = ".";                                                                                              ///< "dot" is a missing string in the VCF spec.
}

inline bool missing (const bool value) { return !value; }                                                                       ///< Returns true if bool is false (missing).
inline bool missing (const float value) { return bcf_float_is_missing(value); }                                                 ///< Returns true if float is missing.
inline bool missing (const int8_t value)  { return value == missing_values::int8; }                                             ///< Returns true if int8_t is missing.
inline bool missing (const int16_t value) { return value == missing_values::int16; }                                            ///< Returns true if int16_t is missing.
inline bool missing (const int32_t value) { return value == missing_values::int32; }                                            ///< Returns true if int32_t is missing.
inline bool missing (const std::string& value) { return value.empty() || value == missing_values::string_dot;}                  ///< Returns true if string is missing.
inline bool missing (const char* value) { return value == missing_values::string_empty || value == missing_values::string_dot;} ///< Returns true if char* is missing.

/**
 * Returns true if value is missing.
 * @tparam MISSING_TYPE any class that implements the missing() as a public member function.
 * @return True if value is missing.
 */
template <class MISSING_TYPE>
inline bool missing(const MISSING_TYPE& value) {
  return value.missing();
}

/**
 * Missing overload for functions that return a vector of values. It only applies if the entire vector is missing.
 * @tparam VALUE any type that can be fit into a container. Any type, really.
 * @param v any vector
 * @return true if the vector is empty (therefore the value that was returned is missing)
 */
template <class VALUE>
inline bool missing(const std::vector<VALUE>& v) { return v.empty(); }

}

#endif // gamgee__missing__guard
