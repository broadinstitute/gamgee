#ifndef gamgee__genotype_field_value_iterator__
#define gamgee__genotype_field_value_iterator__

#include "utils/format_field_type.h"

#include "htslib/vcf.h"

#include <iterator>
#include <memory>

namespace gamgee {

class GenotypeFieldValueIterator : public std::iterator<std::random_access_iterator_tag, int32_t> {
  public:

    GenotypeFieldValueIterator(const std::shared_ptr<bcf1_t>& body, const uint8_t* data_ptr, const uint8_t num_bytes, const utils::FormatFieldType& type) :
      m_body {body},
      m_current_data_ptr {const_cast<uint8_t*>(data_ptr)},
      m_original_data_ptr {const_cast<uint8_t*>(data_ptr)},
      m_num_bytes {num_bytes},
      m_type {type}
    {
    }

    GenotypeFieldValueIterator(const GenotypeFieldValueIterator& other) :
      m_body {other.m_body},
      m_current_data_ptr {other.m_current_data_ptr},
      m_original_data_ptr {other.m_original_data_ptr},
      m_num_bytes {other.m_num_bytes},
      m_type {other.m_type}
    {
    }

    GenotypeFieldValueIterator(GenotypeFieldValueIterator&& other) noexcept :
      m_body {std::move(other.m_body)},
      m_current_data_ptr {other.m_current_data_ptr},
      m_original_data_ptr {other.m_original_data_ptr},
      m_num_bytes {other.m_num_bytes},
      m_type {other.m_type}
    {
    }

    GenotypeFieldValueIterator& operator=(const GenotypeFieldValueIterator& other) {
      if (this == &other)
        return *this;
      m_body = std::move(other.m_body);
      m_current_data_ptr = other.m_current_data_ptr;
      m_original_data_ptr = other.m_original_data_ptr;
      m_num_bytes = other.m_num_bytes;
      m_type = other.m_type;
      return *this;
    }

    GenotypeFieldValueIterator& operator=(GenotypeFieldValueIterator&& other) noexcept {
      if (this == &other)
        return *this;
      m_body = std::move(other.m_body);
      m_current_data_ptr = other.m_current_data_ptr;
      m_original_data_ptr = other.m_original_data_ptr;
      m_num_bytes = other.m_num_bytes;
      m_type = other.m_type;
      return *this;
    }

    GenotypeFieldValueIterator& operator+=(const int n) {
      m_current_data_ptr += n * m_num_bytes;
      return *this;
    }

    GenotypeFieldValueIterator& operator-=(const int n) {
      m_current_data_ptr -= n * m_num_bytes;
      return *this;
    }

    bool operator==(const GenotypeFieldValueIterator& other) {
      return m_body == other.m_body && m_current_data_ptr == other.m_current_data_ptr;
    }

    bool operator!=(const GenotypeFieldValueIterator& other) {
      return m_body != other.m_body || m_current_data_ptr != other.m_current_data_ptr;
    }

    bool operator<(const GenotypeFieldValueIterator& other) {
      return m_body == other.m_body && m_current_data_ptr < other.m_current_data_ptr;
    }

    bool operator>(const GenotypeFieldValueIterator& other) {
      return m_body == other.m_body && m_current_data_ptr > other.m_current_data_ptr;
    }

    bool operator<=(const GenotypeFieldValueIterator& other) {
      return m_body == other.m_body && m_current_data_ptr <= other.m_current_data_ptr;
    }

    bool operator>=(const GenotypeFieldValueIterator& other) {
      return m_body == other.m_body && m_current_data_ptr >= other.m_current_data_ptr;
    }

    int32_t operator*() const noexcept {
      return convert_from_byte_array(m_current_data_ptr, 0);
    }

    int32_t operator++() noexcept {
      m_current_data_ptr += m_num_bytes;
      return convert_from_byte_array(m_current_data_ptr, 0);
    }

    int32_t operator--() {
      m_current_data_ptr -= m_num_bytes;
      return convert_from_byte_array(m_current_data_ptr, 0);
    }

    int32_t operator[](const uint32_t index) const {
      return convert_from_byte_array(m_original_data_ptr, index);
    }

  private:
    std::shared_ptr<bcf1_t> m_body;
    uint8_t* m_current_data_ptr;
    uint8_t* m_original_data_ptr;
    uint8_t m_num_bytes;
    utils::FormatFieldType m_type;

    inline int32_t convert_from_byte_array(const uint8_t* data_ptr, int index) const {
      return utils::convert_data_to_integer(data_ptr, index, m_num_bytes, m_type);
    }
};

};

#endif
