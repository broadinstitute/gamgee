#ifndef gamgee__individual_field_value__guard
#define gamgee__individual_field_value__guard

#include "individual_field_value_iterator.h"

#include "../missing.h"
#include "../utils/hts_memory.h"
#include "../utils/utils.h"
#include "../utils/variant_field_type.h"

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
 * A typical use of the IndividualFieldValue can be exemplified by the
 * phred likelihood accessor in Variant:
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * const auto all_pls = record.integer_individual_field("PL");
 * const auto my_pls = all_pls[2] // take the second sample
 * for_each(my_pls.begin(), my_pls.end(), [](const auto pl) { cout << pl << endl; });
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * IndividualFieldValue objects can also be used in for loops like so: 
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (const auto pl : my_pls) 
 *   cout << pl << endl;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * While the IndividualFieldValue objects are not really intended to be
 * created by the user, they are the returned by the IndividualField iterator for
 * types that don't have a specialized object.  Utilizing them correctly can
 * really simplify your work by leveraging the power of the STL functions.
 *
 * @note all methods are inlined on purpose because they are so simple
 * 
 * @tparam TYPE the object type that holds the values for each sample. For
 * example for GQ it's a uint8_t, some types like GT can have specialized classes
 * like Genotype. For all other types it can be the IndividualFieldValue.
 */
template<class VALUE_TYPE>
class IndividualFieldValue {
 public:

  /**
   * @brief creates a new IndividualFieldValue poiinting to the shared byte array inside the variant object
   * @copydetails IndividualField::IndividualField(const std::shared_ptr<bcf1_t>&, bcf_fmt_t*)
   * @param body the the bcf1_t structure to hold a shared pointer to
   * @param format_ptr the format field pointer inside the body
   * @param data_ptr the location in the specific value inside the format_ptr byte array
   */
  IndividualFieldValue(const std::shared_ptr<bcf1_t>& body, const bcf_fmt_t* const format_ptr, uint8_t* const data_ptr) :
    m_body {body},
    m_format_ptr {format_ptr},
    m_data_ptr {data_ptr},
    m_num_bytes {utils::size_for_type(static_cast<utils::VariantFieldType>(format_ptr->type), format_ptr)}
  {}

  /**
   * @brief copying of the IndividualFieldValue object is not allowed.
   * @copydetails IndividualField::IndividualField(const IndividualField&)
   */
  IndividualFieldValue(const IndividualFieldValue& other) = delete;
  
  /**
   * @brief safely moves the data from one IndividualFieldValue to a new one without making any copies
   * @param other another IndividualFieldValue object
   */
  IndividualFieldValue(IndividualFieldValue&& other) :
    m_body {std::move(other.m_body)},
    m_format_ptr {other.m_format_ptr},
    m_data_ptr {other.m_data_ptr},
    m_num_bytes {other.m_num_bytes}
  {}

  /**
   * @copydoc IndividualField::IndividualFieldValue(const IndividualFieldValue&)
   */
  IndividualFieldValue& operator=(const IndividualFieldValue& other) = delete;

  /**
   * @brief safely moves the data from one IndividualField to the other without making any copies
   * @param other another IndividualFieldValue object
   */
  IndividualFieldValue& operator=(IndividualFieldValue&& other) {
    if (this != &other)
      return *this;
    m_body = std::move(other.m_body);
    m_format_ptr = other.m_format_ptr;
    m_data_ptr = other.m_data_ptr;
    m_num_bytes = other.m_num_bytes;
    return *this;
  }

  /**
   * @brief compares two IndividualFieldValue objects in the following order: memory address, size and values. 
   * @param other something to compare to
   * @return true if the objects are the same (memory address-wise), or contain exactly the same values. Value comparison is dictated by TYPE's operator== implementation
   */
  bool operator==(const IndividualFieldValue& other) const {
    if (this == &other) 
      return true;
    //Use iterators where possible as they take care of field sizes, bcf_*_vector_end
    auto other_iter = other.begin();
    auto other_end = other.end();
    for(const auto& curr_val : *this)
    {
      if(other_iter == other_end)	//different length, this is longer (more valid values) than other
	return false;
      if(!utils::bcf_check_equal_element(curr_val, *other_iter))
	return false;
      ++other_iter;
    }
    //Check if other still has more valid values
    if(other_iter != other_end)
      return false;
    return true;
  }

  /**
   * @brief compares two IndividualFieldValue objects in the following order: memory address, size and values. 
   * @param other something to compare to
   * @return true if the objects are not the same (memory address-wise), or contain different number of values, or the values are not exactly the same. Value comparison is dictated by TYPE's operator== implementation
   */
  bool operator!=(const IndividualFieldValue& other) const {
    return !(*this == other);
  }

  /**
   * @brief random access to the value of a given sample for reading or writing
   * @param index must be between 0 and the number of values per sample for this record 
   * @note implementation guarantees this operation to be O(1)
   * @exception std::out_of_range if sample is out of range
   * @return the value in that index
   */
  VALUE_TYPE operator[](const uint32_t index) const {
    utils::check_max_boundary(index, m_format_ptr->n);
    return convert_from_byte_array(index); 
  }

  /**
   * @brief create a new begin iterator over the values for this sample
   */
  IndividualFieldValueIterator<VALUE_TYPE> begin() const {
    return IndividualFieldValueIterator<VALUE_TYPE>{m_body, m_data_ptr, m_data_ptr + m_format_ptr->size, m_num_bytes, static_cast<utils::VariantFieldType>(m_format_ptr->type)};
  }

  /**
   * @brief create a new end iterator over the values for this sample
   */
  IndividualFieldValueIterator<VALUE_TYPE> end() const {
    return IndividualFieldValueIterator<VALUE_TYPE>{m_body, m_data_ptr + m_format_ptr->size, m_data_ptr + m_format_ptr->size, m_num_bytes, static_cast<utils::VariantFieldType>(m_format_ptr->type)};
  }

  uint32_t size() const { return m_format_ptr->n; } ///< @brief the number of values in this IndividualFieldValue (values per sample)
  VALUE_TYPE front() const { return operator[](0); }                    ///< @brief convenience function to access the first element
  VALUE_TYPE back() const { return operator[](m_format_ptr->n - 1); }   ///< @brief convenience function to access the last element

  /**
   * @brief returns true if all of the values are missing
   */
  bool missing() const {
    // Important: relies on IndividualFieldValueIterator to skip vector end values, if there are any
    for (const auto& value : *this) {
      // use qualifier to avoid recursion
      if (!gamgee::missing(value))
        return false;
    }
    return true;
  }

 private:
  std::shared_ptr<bcf1_t> m_body;
  const bcf_fmt_t* m_format_ptr;
  uint8_t* m_data_ptr;
  uint32_t m_num_bytes;

  VALUE_TYPE convert_from_byte_array(int index) const;
};

/**
 * @brief specialization of the conversion template for int32_t
 */
template<> inline
int32_t IndividualFieldValue<int32_t>::convert_from_byte_array(int index) const {
  return utils::convert_data_to_integer(m_data_ptr, index, m_num_bytes, static_cast<utils::VariantFieldType>(m_format_ptr->type));
}

/**
 * @brief specialization of the conversion template for floats
 */
template<> inline
float IndividualFieldValue<float>::convert_from_byte_array(int index) const {
  return utils::convert_data_to_float(m_data_ptr, index, m_num_bytes, static_cast<utils::VariantFieldType>(m_format_ptr->type));
}

/**
 * @brief specialization of the conversion template for strings
 */
template<> inline
std::string IndividualFieldValue<std::string>::convert_from_byte_array(int index) const {
  return utils::convert_data_to_string(m_data_ptr, index, m_num_bytes, static_cast<utils::VariantFieldType>(m_format_ptr->type));
}

/**
 * @brief  specialization for operator[] for strings
 * String fields have 1 string (at most), it is wrong to have a non-0 index 
 */
template<> inline
std::string IndividualFieldValue<std::string>::operator[](const uint32_t index) const {
  auto is_string_type = utils::is_string_type(m_format_ptr->type);
  auto limit = size();
  auto prefix_msg = "";
  if(is_string_type)
  {
    limit = 1u;
    prefix_msg = "FORMAT fields of type string in VCFs have only 1 element per sample :: ";
  }
  utils::check_max_boundary(index, limit, prefix_msg);
  return convert_from_byte_array(index); 
}

}

#endif // gamgee__individual_field_value__guard
