#ifndef gamgee__variant_field__guard
#define gamgee__variant_field__guard

#include <iostream>

#include "variant_field_iterator.h"
#include "utils/hts_memory.h"
#include "utils/utils.h"

#include "htslib/vcf.h"

#include <memory>
#include <sstream>
#include <stdexcept>

namespace gamgee {

/**
 * @brief A class template to hold the values of a specific Variant's format field for all samples
 *
 * The role of this class is to perform the pointer manipulations behind the scenes that permit the user to 
 * navigate the values of a field without making any copies and benefiting from data locality (all the data
 * is stored contiguously in memory). 
 *
 * For example, the genotype quality field (GQ) will hold one uint8_t value for each sample. The PL field 
 * will hold at least three values per sample. The VariantField will hold all these values per sample in a 
 * random access (and random access iterator compatible) way. The VariantField can be used in any algorithm
 * of the stl that requires random access iterators.
 *
 * A typical use of the VariantField can be examplified by the genotype quality accessor in Variant: 
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * const auto gquals = variant_record.genotype_quals(); 
 * for_each(gquals.begin(), gquals.end(), [](const auto q) { cout << q[0] << endl; });
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * VariantField objects can also be used in for loops like so: 
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (const auto q : variant_record.genotype_quals()) 
 *   cout << q[0] << endl;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * For a field with multiple values, you can also iterate on the resulting VariantFieldValue like so:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * const auto pls = variant_record.phred_likelihoods(); 
 * for(const auto& pl : pls) {
 *   for (const auto x : pl) {
 *     cout << x << endl;  // here x will be hom_ref, het, homvar in each iteration
 *   }
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * While the VariantField objects are not really intended to be created by the user, they are the returned by 
 * many accessors in the Variant API like the genotype_quals() example above. Utilizing them correctly can 
 * really simplify your work by leveraging the power of the STL functions.
 *
 * @note all methods are inlined on purpose because they are so simple
 * 
 * @tparam TYPE the type desired for the given tag. For example for GQ's typically you would request an int32.
 *
 * @warning Currently only int32 and std::string are supported
 * 
 */
template<class TYPE>
class VariantField {
 public:

  /**
   * @brief default constructor of an empty VariantField
   * @warning since private members are immutable, creating an empty VariantField like this means that it will stay empty forever (desired immutability effect)
   * @note empty format fields are created when the field requested is missing in the Variant record
   */
  VariantField() : m_body {nullptr}, m_format_ptr {nullptr} {} ///< for some reason the = default doesn't work here. Need to check why (probably related to the template).

  /**
   * @brief creates a new format field object pointing to the shared byte array inside the Variant object
   * @param body the full Variant object where all the values of the field are stored (shared ownership)
   * @param format_ptr a structure with a pointer to the location in the raw byte array (m_body->indiv) where the format field starts and information about the values (size, type, number,...)
   * @note this implementation doesn't make any copies, simply manages the access to the raw byte array in Variant giving iterator interfaces to the users
   */
  explicit VariantField(const std::shared_ptr<bcf1_t>& body, bcf_fmt_t* format_ptr) :
    m_body {body}, m_format_ptr {format_ptr} 
  {}

  /**
   * @brief copying of the VariantField object is not allowed. Use move constructor instead.
   */
  VariantField(const VariantField& other) = delete;

  /**
   * @brief safely moves the data from one VariantField to a new one without making any copies
   * @param other another VariantField object
   */
  VariantField(VariantField&& other) :
    m_body {std::move(other.m_body)},
    m_format_ptr {other.m_format_ptr}
  {}

  /**
   * @copydoc VariantField::VariantField(const VariantField&)
   */
  VariantField& operator=(const VariantField& other) = delete;

  /**
   * @brief safely moves the data from one VariantField to the other without making any copies
   * @param other another VariantField object
   */
  VariantField& operator=(VariantField&& other) {
    if (this != &other)
      return *this;
    m_body = std::move(other.m_body);
    m_format_ptr = other.m_format_ptr;
    return *this;
  }

  /**
   * @brief random access to the value of a given sample for reading or writing
   * @param sample must be between 0 and the number of samples for this record 
   * @note implementation guarantees this operation to be O(1)
   * @exception std::out_of_range if sample is out of range or entire field is missing and trying to access invalid memory
   * @return the value if it is a basic type (e.g. GQ, GL), or a specific object if it is a complex type (e.g. PL, AD,...)
   */
  TYPE operator[](const uint32_t sample) const {
    if (empty())
      throw std::out_of_range("Tried to index an individual field that is missing with operator[]");
    utils::check_max_boundary(sample, m_body->n_sample);
    return TYPE{m_body, m_format_ptr, m_format_ptr->p + (sample * m_format_ptr->size)}; 
  }

  /**
   * @brief an iterator to the beginning of the object
   */
  VariantFieldIterator<TYPE> begin() const {
    return VariantFieldIterator<TYPE>{m_body, m_format_ptr};
  }

  /**
   * @brief an iterator to the end of the object
   */
  VariantFieldIterator<TYPE> end() const {
    return VariantFieldIterator<TYPE>{m_body, m_format_ptr, true};
  }

  uint32_t size() const { return empty() ? 0 : m_body->n_sample; } ///< @brief the number of values in this VariantField @note this will always be the number of samples in the Variant record, which you can obtain once from it and not repeatedly for every VariantField)
  uint32_t n_samples() const { return size(); }                    ///< @brief just an alias to size() to simplify interfaces
  uint32_t empty() const { return m_body == nullptr; }             ///< @brief checks if the object is empty. @note empty objects are returned when the requested field is missing
  uint32_t missing() const { return empty();}                       ///< @brief checks if the object is empty. @note empty objects are returned when the requested field is missing

 private:
  const std::shared_ptr<bcf1_t> m_body; ///< shared ownership of the Variant record memory so it stays alive while this object is in scope
  const bcf_fmt_t* const m_format_ptr;  ///< pointer to m_body structure where the data for this particular type is located.
};


}

#endif // gamgee__variant_field__guard
