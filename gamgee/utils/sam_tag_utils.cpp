#include "sam_tag_utils.h"

#include <cstring>

using namespace std;

namespace gamgee {
namespace utils {


uint8_t* skip_numeric_tag_value(uint8_t* buffer, uint32_t size) {
  return buffer + size;
}

uint8_t* skip_string_tag_value(uint8_t* buffer) {
  while (*buffer) {
    ++buffer;
  }
  return ++buffer;
}

uint8_t* skip_array_tag_value(uint8_t* buffer, uint32_t size) {
  int32_t elements;
  memcpy(&elements, ++buffer, 4);
  buffer += 4;
  return buffer + size * elements;
}

std::unordered_map<std::string, uint8_t*> parse_encoded_tags(uint8_t* buffer, int32_t len) {
  auto result = unordered_map<std::string, uint8_t*>{};
  auto buffer_end = buffer + len;
  while (buffer < buffer_end) {
    string tag_name(reinterpret_cast<char *>(buffer), 2);
    buffer += 2;
    auto t = reinterpret_cast<unsigned char>(*buffer++);
    switch (t) {
      case CHAR_SAM_TAG_TYPE_CODE:
        result[tag_name] = buffer-1;
        buffer = skip_numeric_tag_value(buffer, 1);
        break;
      case INT8_SAM_TAG_TYPE_CODE: case UINT8_SAM_TAG_TYPE_CODE:
        result[tag_name] = buffer-1;
        buffer = skip_numeric_tag_value(buffer, 1);
        break;
      case INT16_SAM_TAG_TYPE_CODE: case UINT16_SAM_TAG_TYPE_CODE:
        result[tag_name] = buffer-1;
        buffer = skip_numeric_tag_value(buffer, 2);
        break;
      case INT32_SAM_TAG_TYPE_CODE: case UINT32_SAM_TAG_TYPE_CODE:
        result[tag_name] = buffer-1;
        buffer = skip_numeric_tag_value(buffer, 4);
        break;
      case FLOAT_SAM_TAG_TYPE_CODE:
        result[tag_name] = buffer-1;
        buffer = skip_numeric_tag_value(buffer, 4);
        break;
      case DOUBLE_SAM_TAG_TYPE_CODE:
        result[tag_name] = buffer-1;
        buffer = skip_numeric_tag_value(buffer, 8);
        break;
      case STRING_SAM_TAG_TYPE_CODE:
        result[tag_name] = buffer-1;
        buffer = skip_string_tag_value(buffer);
        break;
      case BYTE_ARRAY_SAM_TAG_TYPE_CODE:
        result[tag_name] = buffer-1;
        buffer = skip_string_tag_value(buffer);
        break;
      case NUMERIC_ARRAY_SAM_TAG_TYPE_CODE:
        result[tag_name] = buffer-1;
        auto array_t = reinterpret_cast<unsigned char>(*buffer);
        switch (array_t) {
          case INT8_SAM_TAG_TYPE_CODE: case UINT8_SAM_TAG_TYPE_CODE:
            buffer = skip_array_tag_value(buffer, 1);
            break;
          case INT16_SAM_TAG_TYPE_CODE: case UINT16_SAM_TAG_TYPE_CODE:
            buffer = skip_array_tag_value(buffer, 2);
            break;
          case INT32_SAM_TAG_TYPE_CODE: case UINT32_SAM_TAG_TYPE_CODE:
            buffer = skip_array_tag_value(buffer, 4);
            break;
          case FLOAT_SAM_TAG_TYPE_CODE:
            buffer = skip_array_tag_value(buffer, 4);
            break;
        }
        break;
   }
 }
  return result;
}

SamTagType numeric_array_tag_type(const char& type_char) {
  switch (type_char) {
    case INT8_SAM_TAG_TYPE_CODE:   return SamTagType::INTEGER8ARRAY;
    case UINT8_SAM_TAG_TYPE_CODE:  return SamTagType::UINTEGER8ARRAY;
    case INT16_SAM_TAG_TYPE_CODE:  return SamTagType::INTEGER16ARRAY;
    case UINT16_SAM_TAG_TYPE_CODE: return SamTagType::UINTEGER16ARRAY;
    case INT32_SAM_TAG_TYPE_CODE:  return SamTagType::INTEGER32ARRAY;
    case UINT32_SAM_TAG_TYPE_CODE: return SamTagType::UINTEGER32ARRAY;
    case FLOAT_SAM_TAG_TYPE_CODE:  return SamTagType::FLOATARRAY;
    default: throw invalid_argument("Unsupported sam numeric array tag type: " +
                                    std::to_string(type_char));
  }
}

char numeric_array_tag_type(const SamTagType& type) {
  switch (type) {
    case SamTagType::INTEGER8ARRAY:   return INT8_SAM_TAG_TYPE_CODE;
    case SamTagType::UINTEGER8ARRAY:  return UINT8_SAM_TAG_TYPE_CODE;
    case SamTagType::INTEGER16ARRAY:  return INT16_SAM_TAG_TYPE_CODE;
    case SamTagType::UINTEGER16ARRAY: return UINT16_SAM_TAG_TYPE_CODE;
    case SamTagType::INTEGER32ARRAY:  return INT32_SAM_TAG_TYPE_CODE;
    case SamTagType::UINTEGER32ARRAY: return UINT32_SAM_TAG_TYPE_CODE;
    case SamTagType::FLOATARRAY:      return FLOAT_SAM_TAG_TYPE_CODE;
    default: throw invalid_argument("Unsupported sam numeric array tag type.");
  }
}

}
}
