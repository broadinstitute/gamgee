#ifndef gamgee__sam_tag__guard
#define gamgee__sam_tag__guard

#include <vector>
#include <string>

namespace gamgee {

/**
 * @brief comprehensive list of valid tag types
 *
 * @note CHARARRAY, INTEGERARRAY, and DOUBLEARRAY enums are for backward compatibility.
 * Sam::all_tags() (note: deprecated) returns aux array types ('B') as
 *   CHARARRAY (corresponding to int8, uint8),
 *   INTEGERARRAY (corresponding to int16, uint16, int32, uint32),
 *   and DOUBLEARRAY (corresponds to float).
 * Sam::all_tag_types() use more detailed aux array type enums:
 *   [U]INTEGER8ARRAY, [U]INTEGER16ARRAY, [U]INTEGER32ARRAY, and FLOATARRAY.
 */
enum class SamTagType {
  NONE = -1,
  INTEGER,
  FLOAT,
  DOUBLE,
  CHAR,
  STRING,
  BYTEARRAY,
  CHARARRAY,       // Note: DEPRECATED. only used in Sam::all_tags()
  INTEGERARRAY,    // Note: DEPRECATED. only used in Sam::all_tags()
  INTEGER8ARRAY,
  UINTEGER8ARRAY,
  INTEGER16ARRAY,
  UINTEGER16ARRAY,
  INTEGER32ARRAY,
  UINTEGER32ARRAY,
  FLOATARRAY,
  DOUBLEARRAY     // Note: DEPRECATED. only used in Sam::all_tags()
};

constexpr char CHAR_SAM_TAG_TYPE_CODE = 'A';
constexpr char INTEGER_SAM_TAG_TYPE_CODE = 'i';
constexpr char FLOAT_SAM_TAG_TYPE_CODE = 'f';
constexpr char DOUBLE_SAM_TAG_TYPE_CODE = 'd';
constexpr char STRING_SAM_TAG_TYPE_CODE = 'Z';
constexpr char BYTE_ARRAY_SAM_TAG_TYPE_CODE = 'H';
constexpr char NUMERIC_ARRAY_SAM_TAG_TYPE_CODE = 'B';
// Numeric array types
constexpr char INT8_SAM_TAG_TYPE_CODE = 'c';
constexpr char UINT8_SAM_TAG_TYPE_CODE = 'C';
constexpr char INT16_SAM_TAG_TYPE_CODE = 's';
constexpr char UINT16_SAM_TAG_TYPE_CODE = 'S';
constexpr char INT32_SAM_TAG_TYPE_CODE = 'i';
constexpr char UINT32_SAM_TAG_TYPE_CODE = 'I';

constexpr int32_t SAM_TAG_NAME_LENGTH = 2;
constexpr int32_t SAM_TAG_TYPE_LENGTH = 1;

/**
 * @brief class to represent a Sam TAG:TYPE:VALUE entry
 */
template<class TAG_TYPE>
class SamTag {
 public:
  explicit SamTag(const std::string& name, const TAG_TYPE& value, const bool missing = false) :
    m_name { name },
    m_missing { missing }
  { m_value = value; }

  explicit SamTag(const std::string& name, TAG_TYPE&& value, const bool missing = false) :
    m_name { name },
    m_value { std::move(value) },
    m_missing { missing }
  {}

  SamTag(const SamTag& other) = default;
  SamTag(SamTag&& other) = default;
  SamTag& operator=(const SamTag& other) = default;
  SamTag& operator=(SamTag&& other) = default;
  ~SamTag() = default;

  std::string name() const { return m_name; }
  TAG_TYPE value() const { return m_value; }
  bool missing() const { return m_missing; }

 private:
  std::string m_name;
  TAG_TYPE m_value;
  bool m_missing;
};

/**
 * @brief class to represent sam numeric array tags.
 */
class SamNumericArrayTag {
 public:
  SamNumericArrayTag() :
    m_type { SamTagType::NONE },
    m_integer_values {},
    m_float_values {}
  {};
  SamNumericArrayTag(const SamTagType& type, const std::vector<int64_t>& values);
  SamNumericArrayTag(const SamTagType& type, const std::vector<float>& values);
  SamNumericArrayTag(const SamNumericArrayTag& other) = default;
  ~SamNumericArrayTag() = default;

  bool operator==(const SamNumericArrayTag&  numeric_array_tag) const;

  SamTagType type() const { return m_type; }
  std::vector<int64_t> int_values() const { return m_integer_values; }
  std::vector<float> float_values() const { return m_float_values; }
  int64_t int_value(const uint32_t index) const { return m_integer_values[index]; }
  float float_value(const uint32_t index) const { return m_float_values[index]; }
  uint32_t size() const {
    if ( m_type == SamTagType::FLOATARRAY ) return m_float_values.size();
    else return m_integer_values.size();
  }

 private:
  // One of [U]INTEGER8ARRAY, [U]INTEGER16ARRAY, [U]INTEGER32ARRAY, and FLOATARRAY ).
  SamTagType m_type;
  // To hold values for array type [U]INTEGER8ARRAY, or [U]INTEGER16ARRAY, or [U]INTEGER32ARRAY
  std::vector<int64_t> m_integer_values;
  // To hold values for array type FLOATARRAY
  std::vector<float> m_float_values;
};

}


#endif // gamgee__sam_tag__guard
