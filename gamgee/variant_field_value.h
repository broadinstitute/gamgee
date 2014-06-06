#ifndef __gamgee__variant_field_value__
#define __gamgee__variant_field_value__

#include "variant_field_value_iterator.h"
#include "utils/hts_memory.h"
#include "utils/utils.h"
#include "utils/format_field_type.h"

#include "htslib/vcf.h"

#include <memory>

namespace gamgee {

/**
 * @brief A class template to hold the values a format field for a particular sample
 *
 * The role of this class is to perform the pointer manipulations behind the
 * scenes that permit the user to navigate the values of a field in a sample
 * without making any copies and benefiting from data locality (all the data is
 * stored contiguously in memory). 
 *
 * A typical use of the VariantFieldValue can be examplified by the
 * genotype qualitiy accessor in Variant: 
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * const auto all_pls = record.phred_likelihoods();
 * const auto my_pls = all_pls[2] // take the second sample
 * for_each(my_pls.begin(), my_pls.end(), [](const auto pl) { cout << pl << endl; });
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * VariantFieldValue objects can also be used in for loops like so: 
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (const auto pl : my_pls) 
 *   cout << pl << endl;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * While the VariantFieldValue objects are not really intended to be
 * created by the user, they are the returned by the VariantField iterator for
 * types that don't have a specialized object.  Utilizing them correctly can
 * really simplify your work by leveraging the power of the STL functions.
 *
 * @note all methods are inlined on purpose because they are so simple
 * 
 * @tparam TYPE the object type that holds the values for each sample. For
 * example for GQ it's a uint8_t, some types like GT and PL can have
 * specialized classes like Genotype and PhredLikelihoods. For all other types
 * it can be the VariantFieldValue.
 */
template<class VALUE_TYPE>
class VariantFieldValue {
 public:

  /**
   * @brief creates a new VariantFieldValue poiinting to the shared byte array inside the variant object
   * @copydetails VariantField::VariantField(const std::shared_ptr<bcf1_t>&, bcf_fmt_t*)
   * @param body the the bcf1_t structure to hold a shared pointer to
   * @param format_ptr the format field pointer inside the body
   * @param data_ptr the location in the specific value inside the format_ptr byte array
   */
  VariantFieldValue(const std::shared_ptr<bcf1_t>& body, const bcf_fmt_t* const format_ptr, uint8_t* const data_ptr) :
    m_body {body},
    m_format_ptr {format_ptr},
    m_data_ptr {data_ptr},
    m_num_bytes {utils::size_for_type(static_cast<utils::FormatFieldType>(format_ptr->type))}
  {}

  /**
   * @brief copying of the VariantFieldValue object is not allowed.
   * @copydetails VariantField::VariantField(const VariantField&)
   */
  VariantFieldValue(const VariantFieldValue& other) = delete;
  
  /**
   * @brief safely moves the data from one VariantFieldValue to a new one without making any copies
   * @param other another VariantFieldValue object
   */
  VariantFieldValue(VariantFieldValue&& other) :
    m_body {std::move(other.m_body)},
    m_format_ptr {other.m_format_ptr},
    m_data_ptr {other.m_data_ptr},
    m_num_bytes {other.m_num_bytes}
  {}

  /**
   * @copydoc VariantField::VariantFieldValue(const VariantFieldValue&)
   */
  VariantFieldValue& operator=(const VariantFieldValue& other) = delete;

  /**
   * @brief safely moves the data from one VariantField to the other without making any copies
   * @param other another VariantFieldValue object
   */
  VariantFieldValue& operator=(VariantFieldValue&& other) {
    if (this != &other)
      return *this;
    m_body = std::move(other.m_body);
    m_format_ptr = other.m_format_ptr;
    m_data_ptr = other.m_data_ptr;
    m_num_bytes = other.m_num_bytes;
    return *this;
  }

  /**
   * @brief random access to the value of a given sample for reading or writing
   * @param index must be between 0 and the number of values per sample for this record 
   * @note implementation guarantees this operation to be O(1)
   * @exception std::out_of_range if sample is out of range
   * @return the value in that index
   */
  VALUE_TYPE operator[](const uint32_t index) const {
    utils::check_boundaries(index, m_format_ptr->n - 1);
    return convert_from_byte_array(index); 
  }

  /**
   * @brief create a new begin iterator over the values for this sample
   */
  VariantFieldValueIterator<VALUE_TYPE> begin() const {
    return VariantFieldValueIterator<VALUE_TYPE>{m_body, m_data_ptr, m_num_bytes, static_cast<utils::FormatFieldType>(m_format_ptr->type)};
  }

  /**
   * @brief create a new end iterator over the values for this sample
   */
  VariantFieldValueIterator<VALUE_TYPE> end() const {
    return VariantFieldValueIterator<VALUE_TYPE>{m_body, m_data_ptr + m_format_ptr->size, m_num_bytes, static_cast<utils::FormatFieldType>(m_format_ptr->type)};
  }

 private:
  const std::shared_ptr<bcf1_t> m_body;
  const bcf_fmt_t* const m_format_ptr;
  uint8_t* const m_data_ptr;
  const uint8_t m_num_bytes;

  VALUE_TYPE convert_from_byte_array(int index) const;
};

/**
 * @brief specialization of the conversion template for int32_t
 */
template<> inline
int32_t VariantFieldValue<int32_t>::convert_from_byte_array(int index) const {
  return utils::convert_data_to_integer(m_data_ptr, index, m_num_bytes, static_cast<utils::FormatFieldType>(m_format_ptr->type));
}

/**
 * @brief specialization of the conversion template for floats
 */
template<> inline
float VariantFieldValue<float>::convert_from_byte_array(int index) const {
  return utils::convert_data_to_float(m_data_ptr, index, m_num_bytes, static_cast<utils::FormatFieldType>(m_format_ptr->type));
}

/**
 * @brief specialization of the conversion template for strings
 */
template<> inline
std::string VariantFieldValue<std::string>::convert_from_byte_array(int index) const {
  return utils::convert_data_to_string(m_data_ptr, index, m_num_bytes, static_cast<utils::FormatFieldType>(m_format_ptr->type));
}


}

#endif
