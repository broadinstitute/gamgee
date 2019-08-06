#ifndef gamgee__sam_tag_utils__guard
#define gamgee__sam_tag_utils__guard

#include "../sam/sam_tag.h"

#include <string>
#include <unordered_map>

namespace gamgee {

namespace utils {

/**
 * @brief skips the bytes representing the tag's value for numeric tags.
 * @param buffer pointer to the first byte in the tag value
 */
uint8_t* skip_numeric_tag_value(uint8_t* buffer, uint32_t size);

/**
 * @brief skips the bytes representing the tag's value for string tags.
 * @param buffer pointer to the first byte in the tag value string
 */
uint8_t* skip_string_tag_value(uint8_t* buffer);

/**
 * @brief skips the bytes representing the tag's value for array tags.
 * @param buffer pointer to the first byte in the tag value array
 */
uint8_t* skip_array_tag_value(uint8_t* buffer, uint32_t size);

/**
 * @brief parse input encoded tag record and retrieve all tag names and pointer to their type:value
 */
std::unordered_map<std::string, uint8_t*> parse_encoded_tags(uint8_t* buffer, int32_t len);

/**
 * @brief return the type of sam numeric array tag.
 */
SamTagType numeric_array_tag_type(const char& type_char);

/**
 * @brief return the character representing the type of sam numeric array tag.
 */
char numeric_array_tag_type(const SamTagType& type);

} // end utils namespace
} // end gamgee namespace

#endif // gamgee__sam_tag_utils__guard
