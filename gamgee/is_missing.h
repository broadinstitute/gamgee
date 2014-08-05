#ifndef gamgee__is_missing__guard
#define gamgee__is_missing__guard

#include "sam_tag.h"
#include "htslib/vcf.h"

#include <string>
#include <cmath>

namespace gamgee {

namespace missing_values {
  constexpr auto int8 = bcf_int8_missing;                                                                                          ///< missing value for an int8
  constexpr auto int16 = bcf_int16_missing;                                                                                        ///< missing value for an int16
  constexpr auto int32 = bcf_int32_missing;                                                                                        ///< missing value for an int32
  constexpr auto string_empty = "";                                                                                                ///< empty string is a missing string
  constexpr auto string_dot = ".";                                                                                                 ///< "dot" is a missing string in the VCF spec.
}

inline bool is_missing (const float value) { return isnan(value); }                                                                ///< @brief checks whether or not a float value returned by a Variant field is missing
inline bool is_missing (const int8_t value)  { return value == missing_values::int8; }                                             ///< @brief Returns true if value is missing.
inline bool is_missing (const int16_t value) { return value == missing_values::int16; }                                            ///< @brief Returns true if value is missing.
inline bool is_missing (const int32_t value) { return value == missing_values::int32; }                                            ///< @brief Returns true if value is missing.
inline bool is_missing (const std::string& value) { return value.empty() || value == missing_values::string_dot;}                  ///< @brief Returns true if value is missing.
inline bool is_missing (const char* value) { return value == missing_values::string_empty || value == missing_values::string_dot;} ///< @brief Returns true if value is missing.

/**
 * @brief Returns true if value is missing.
 * @tparam MISSING_TYPE any class that implements the is_missing() as a public member function.
 * @return True if value is missing.
 */
template <class MISSING_TYPE>
inline bool is_missing(const MISSING_TYPE& value) {
  return value.is_missing();
}

}

#endif // gamgee__is_missing__guard
