
#include "sam_tag.h"

#include <stdexcept>

using namespace std;

namespace gamgee {

/**
 * @brief creates a SamNumericArrayTag object of integer type.
*/
SamNumericArrayTag::SamNumericArrayTag(const SamTagType& type,
                                       const std::vector<int64_t>& values) :
  m_type { type },
  m_integer_values { values },
  m_float_values {}
{
  if ( (type >= SamTagType::UINTEGER32ARRAY) || (type <= SamTagType::INTEGER8ARRAY) )
      invalid_argument("Numeric array type is not integer.");
}

/**
 * @brief creates a SamNumericArrayTag object of float type.
*/
SamNumericArrayTag::SamNumericArrayTag(const SamTagType& type,
                                       const std::vector<float>& values) :
  m_type { type },
  m_integer_values {},
  m_float_values { values }
{
  if ( type != SamTagType::FLOATARRAY )
      invalid_argument("Numeric array type is not float.");
}

/**
 * @brief defines "==" operator for SamNumericArrayTag object.
*/
bool SamNumericArrayTag::operator==(const SamNumericArrayTag&  numeric_array_tag) const {
  bool result = (numeric_array_tag.type() == m_type) &&
                (this->size() == numeric_array_tag.size());
  if (!result) return result;

  if (m_type == SamTagType::FLOATARRAY) {
    for (uint32_t i = 0; i < this->size(); ++i) {
      result = result &&
              ( this->float_value(i) == numeric_array_tag.float_value(i) );
      if ( !result ) return result;
    }
  } else {  // {U}INTERGER{8,16,32}ARRAY type
    for (uint32_t i = 0; i < this->size(); ++i) {
      result = result &&
              ( this->int_value(i) == numeric_array_tag.int_value(i) );
      if ( !result ) return result;
    }
  }

  return result;
}

}
