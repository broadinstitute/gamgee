#ifndef gamgee__shared_field__guard
#define gamgee__shared_field__guard

#include "shared_field_iterator.h"

#include "../utils/utils.h"

#include "htslib/vcf.h"

#include <memory>
#include <sstream>
#include <stdexcept>

namespace gamgee {

/**
 * @brief A class template to hold the values of a specific Variant's shared field 
 *
 * The role of this class is to perform the pointer manipulations behind the scenes that permit the user to 
 * navigate the values of a field without making any copies and benefiting from data locality (all the data
 * is stored contiguously in memory). 
 *
 * For example, the allele number shared field (AN) will hold one integer value. A ranksum test field for example
 * could hold a vector of floats. The SharedField object will hold all these values for a shared field in a 
 * random access (and random access iterator compatible) way. The SharedField can be used in any algorithm
 * of the STL that requires random access iterators.
 *
 * A typical use of the SharedField can be exemplified by the following:
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * const auto an = variant_record.shared_integer_field("AN"); 
 * for_each(an.begin(), an.end(), [](const auto x) { cout << x[0] << endl; });
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * SharedField objects can also be used in for loops like so: 
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (const auto q : variant_record.shared_string_field("CULPRIT")) 
 *   cout << q << endl;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Or directly like so: 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * const auto ab = variant_record.shared_float_field("AB");
 * cout << ab[0] << endl;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *
 * While the SharedField objects are not really intended to be created by the user, they are returned by 
 * many accessors in the Variant API like the Variant::shared_integer_field() example above. Utilizing them correctly can 
 * really simplify your work by leveraging the power of the STL functions.
 *
 * @note all methods are inlined on purpose because they are so simple
 * @tparam TYPE the output type desired for the given tag. For example for AN's typically you would request an int32_t. Conversions are allowed if possible.
 */
template<class TYPE>
class SharedField {
 public:

  /**
   * @brief default constructor of an empty SharedField
   * @warning since private members are immutable, creating an empty SharedField like this means that it will stay empty forever (desired immutability effect)
   * @note empty shared fields are created when the field requested is missing in the Variant record
   */
  SharedField() : m_body {nullptr}, m_info_ptr {nullptr}, m_bytes_per_value {0} {} 

  /**
   * @brief creates a new SharedField pointing to the shared byte array inside the variant object
   * @param body the the bcf1_t structure to hold a shared pointer to
   * @param info_ptr the info field pointer inside the body
   */
  explicit SharedField(const std::shared_ptr<bcf1_t>& body, const bcf_info_t* const info_ptr) :
    m_body {body},
    m_info_ptr {info_ptr},
    m_bytes_per_value {utils::size_for_type(static_cast<utils::VariantFieldType>(info_ptr->type), info_ptr)}
  {}

  SharedField(const SharedField& other) = delete;            ///< @brief copying of the SharedField object is not allowed. Use move constructor instead.
  SharedField& operator=(const SharedField& other) = delete; ///< @copydoc SharedField::SharedField(const SharedField&)
  SharedField(SharedField&& other) = default;                ///< @brief safely moves the data from one SharedField to a new one without making any copies
  SharedField& operator=(SharedField&& other) = default;     ///< @brief safely moves the data from one SharedField to the other without making any copies

  /**
   * @brief compares two SharedField objects in the following order: memory address, size and values. 
   * @param other something to compare to
   * @return true if the objects are the same (memory address-wise), or contain exactly the same values. Value comparison is dictated by TYPE's operator== implementation
   */
  bool operator==(const SharedField& other) const {
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
   * @brief compares two SharedField objects in the following order: memory address, size and values. 
   * @param other something to compare to
   * @return true if the objects are not the same (memory address-wise), or contain different number of values, or the values are not exactly the same. Value comparison is dictated by TYPE's operator== implementation
   */
  bool operator!=(const SharedField& other) const {
    return !(*this == other);
  }


  /**
   * @brief random access to a given value for reading or writing
   * @param index must be between 0 and the number of values for this record 
   * @note implementation guarantees this operation to be O(1)
   * @exception std::out_of_range if index is out of range or entire field is missing and trying to access invalid memory
   * @return the value in that index
   */
  TYPE operator[](const uint32_t index) const {
    if (empty())
      throw std::out_of_range("Tried to index a shared field that is missing with operator[]");
    utils::check_max_boundary(index, m_info_ptr->len);
    return convert_from_byte_array(index); 
  }

  /** @brief create a new iterator pointing to the begining of the values of this field */
  SharedFieldIterator<TYPE> begin() const {
    if (empty()) return SharedFieldIterator<TYPE>();
    return SharedFieldIterator<TYPE>{m_body, m_info_ptr->vptr, m_info_ptr->vptr + size() * m_bytes_per_value, m_bytes_per_value, static_cast<utils::VariantFieldType>(m_info_ptr->type)};
  }

  /** @brief create a new iterator pointing to the end of the values of this field */
  SharedFieldIterator<TYPE> end() const {
    if (empty()) return SharedFieldIterator<TYPE>();
    return SharedFieldIterator<TYPE>{m_body, m_info_ptr->vptr + size() * m_bytes_per_value, m_bytes_per_value, static_cast<utils::VariantFieldType>(m_info_ptr->type)};
  }

  uint32_t size() const { return m_info_ptr->len; } ///< @brief the number of values in this SharedField (values per sample)
  bool empty() const { return m_body == nullptr; }  ///< @brief checks if the object is empty. @note empty objects are returned when the requested field is missing
  bool missing() const { return empty(); }          ///< simple overload of empty to work with gamgee's missing API
  TYPE front() const { return operator[](0); }                    ///< @brief convenience function to access the first element
  TYPE back() const { return operator[](m_info_ptr->len - 1); }   ///< @brief convenience function to access the last element

 private:
  std::shared_ptr<bcf1_t> m_body;
  const bcf_info_t* m_info_ptr;
  uint32_t m_bytes_per_value;

  TYPE convert_from_byte_array(int index) const;
};

template<> inline
uint32_t SharedField<std::string>::size() const { return empty() ? 0 : 1; } ///< string info fields cannot take multiple values (by spec)

/**
 * @brief specialization of the conversion template for int32_t
 */
template<> inline
int32_t SharedField<int32_t>::convert_from_byte_array(int index) const {
  return utils::convert_data_to_integer(m_info_ptr->vptr, index, m_bytes_per_value, static_cast<utils::VariantFieldType>(m_info_ptr->type));
}

/**
 * @brief specialization of the conversion template for floats
 */
template<> inline
float SharedField<float>::convert_from_byte_array(int index) const {
  return utils::convert_data_to_float(m_info_ptr->vptr, index, m_bytes_per_value, static_cast<utils::VariantFieldType>(m_info_ptr->type));
}

/**
 * @brief specialization of the conversion template for strings
 */
template<> inline
std::string SharedField<std::string>::convert_from_byte_array(int index) const {
  return utils::convert_data_to_string(m_info_ptr->vptr, index, m_bytes_per_value, static_cast<utils::VariantFieldType>(m_info_ptr->type));
}

/**
 * @brief  specialization for operator[] for strings
 * String fields have 1 string (at most), it is wrong to have a non-0 index 
 */
template<> inline
std::string SharedField<std::string>::operator[](const uint32_t index) const {
  if (empty())
    throw std::out_of_range("Tried to index a shared field that is missing with operator[]");
  auto is_string_type = utils::is_string_type(m_info_ptr->type);
  auto limit = m_info_ptr->len; //cannot use size here as size<string> is specialized to return 1
  auto prefix_msg = "";
  if(is_string_type)
  {
    limit = 1u;
    prefix_msg = "INFO fields of type string in VCFs have only 1 element per sample :: ";
  }
  utils::check_max_boundary(index, limit, prefix_msg);
  return convert_from_byte_array(index); 
}

} // end of namespace gamgee

#endif // gamgee__shared_field__guard
