#include "sam_builder_data_field.h"

#include <cstring>

using namespace std;

namespace gamgee {

/**
 * @brief initialize a SamBuilderDataField to an empty value
 */
SamBuilderDataField::SamBuilderDataField() :
  m_data {},
  m_num_bytes { 0 },  // Note: default no-arg constructor would NOT zero out the POD members
  m_num_elements { 0 }
{}

/**
 * @brief initialize a SamBuilderDataField by copying data from a raw pointer
 *
 * @note takes no ownership of copy_source
 */
SamBuilderDataField::SamBuilderDataField(const void* copy_source, const uint32_t bytes_to_copy, const uint32_t num_elements) :
  m_data { new uint8_t[bytes_to_copy] },
  m_num_bytes { bytes_to_copy },
  m_num_elements { num_elements }
{
  memcpy(m_data.get(), copy_source, bytes_to_copy);
}

/**
 * @brief initialize a SamBuilderDataField by moving an existing unique_ptr into it (without copying the existing data)
 *
 * @note takes ownership of the memory managed by move_source
 */
SamBuilderDataField::SamBuilderDataField(std::unique_ptr<uint8_t[]>&& move_source, const uint32_t source_bytes, const uint32_t num_elements) :
  m_data { move(move_source) },
  m_num_bytes { source_bytes },
  m_num_elements { num_elements }
{}

/**
 * @brief initialize a SamBuilderDataField via move from an existing field
 */
SamBuilderDataField::SamBuilderDataField(SamBuilderDataField&& other) :
  m_data { move(other.m_data) },
  m_num_bytes { other.m_num_bytes },
  m_num_elements { other.m_num_elements }
{}

/**
 * @brief move an existing SamBuildDataField into this one
 */
SamBuilderDataField& SamBuilderDataField::operator=(SamBuilderDataField&& other) {
  if ( &other == this )
    return *this;

  m_data = move(other.m_data);
  m_num_bytes = other.m_num_bytes;
  m_num_elements = other.m_num_elements;
  return *this;
}

/**
 * @brief update the field by copying data from a raw pointer (takes no ownership of copy_source)
 *
 * @note any previous value of m_data is destroyed via the unique_ptr assignment
 */
void SamBuilderDataField::update(const void* copy_source, const uint32_t bytes_to_copy, const uint32_t num_elements) {
  if (num_elements == 0) {
    if (bytes_to_copy != 0)
      throw invalid_argument(string("Non-zero bytes_to_copy for empty cigar."));
    m_data = unique_ptr<uint8_t[]>{ new uint8_t[bytes_to_copy] };
  } else {
    m_data = unique_ptr<uint8_t[]>{ new uint8_t[bytes_to_copy] };
    memcpy(m_data.get(), copy_source, bytes_to_copy);
  }
  m_num_bytes = bytes_to_copy;
  m_num_elements = num_elements;
}

/**
 * @brief update the field by moving an existing unique_ptr into it and taking ownership (without copying the existing data)
 *
 * @note any previous value of m_data is destroyed via the unique_ptr assignment
 */
void SamBuilderDataField::update(std::unique_ptr<uint8_t[]>&& move_source, const uint32_t source_bytes, const uint32_t num_elements) {
  m_data = move(move_source);
  m_num_bytes = source_bytes;
  m_num_elements = num_elements;
}

/**
 * @brief copy this field's byte array into an arbitrary location
 *
 * @return pointer to the byte just AFTER the end of the copied data
 */
uint8_t* SamBuilderDataField::copy_into(uint8_t* destination) const {
  memcpy(destination, m_data.get(), m_num_bytes);
  return destination + m_num_bytes;
}

}
