#ifndef __gamgee__sam_builder_data_field__
#define __gamgee__sam_builder_data_field__

#include <memory>

namespace gamgee {

/**
 * @brief class to hold encoded byte arrays for individual data fields (cigar, bases, etc.) during building of a Sam
 *
 * Fields can be created/updated either by copying data from raw pointers, or by moving managed pointers
 * into them without copying. Eg.,
 *
 * auto field = SamBuilderDataField{raw_pointer, num_bytes, num_elements};      // does a copy
 * auto field = SamBuilderDataField{move(unique_ptr), num_bytes, num_elements}; // no copy
 *
 * After construction, a field's value can be altered via the update() functions
 */
class SamBuilderDataField {
 public:
  SamBuilderDataField();                                                                                                             ///<  @brief initialize a SamBuilderDataField to an empty value
  explicit SamBuilderDataField(const void* copy_source, const uint32_t bytes_to_copy, const uint32_t num_elements);                  ///<  @brief initialize a SamBuilderDataField by copying data from a raw pointer (takes no ownership of copy_source)
  explicit SamBuilderDataField(std::unique_ptr<uint8_t[]>&& move_source, const uint32_t source_bytes, const uint32_t num_elements);  ///<  @brief initialize a SamBuilderDataField by moving an existing unique_ptr into it and taking ownership (without copying the existing data)
  SamBuilderDataField(SamBuilderDataField&& other);
  SamBuilderDataField& operator=(SamBuilderDataField&& other);

  // Fields cannot be copied (only moved), and use default destruction
  SamBuilderDataField(const SamBuilderDataField& other) = delete;
  SamBuilderDataField& operator=(const SamBuilderDataField& other) = delete;
  ~SamBuilderDataField() = default;

  const uint8_t* raw_data_ptr() const { return m_data.get(); }  ///< gets a raw pointer to the data buffer
  uint32_t num_bytes() const { return m_num_bytes; }            ///< number of bytes in the data buffer
  uint32_t num_elements() const { return m_num_elements; }      ///< number of elements (cigar operations, bases, etc.) in the data buffer
  bool is_empty() const { return m_num_bytes == 0; }            ///< does this field have any data?

  void update(const void* copy_source, const uint32_t bytes_to_copy, const uint32_t num_elements);                   ///<  @brief update the field by copying data from a raw pointer (takes no ownership of copy_source)
  void update(std::unique_ptr<uint8_t[]>&& move_source, const uint32_t source_bytes, const uint32_t num_elements);   ///<  @brief update the field by moving an existing unique_ptr into it and taking ownership (without copying the existing data)
  uint8_t* copy_into(uint8_t* destination) const;    ///< @brief copy this field's byte array into an arbitrary location

 private:
  std::unique_ptr<uint8_t[]> m_data;  ///< buffer containing encoded data for the field, managed exclusively by us
  uint32_t m_num_bytes;               ///< number of bytes in m_data
  uint32_t m_num_elements;            ///< number of elements (cigar operations, bases, etc.) in m_data
};

}

#endif /* __gamgee__sam_builder_data_field__ */
