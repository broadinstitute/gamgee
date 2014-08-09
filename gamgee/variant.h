#ifndef gamgee__variant__guard
#define gamgee__variant__guard

#include "variant_header.h"
#include "individual_field.h"
#include "individual_field_value.h"
#include "shared_field.h"
#include "variant_filters.h"
#include "boost/dynamic_bitset.hpp"
#include "genotype.h"

#include "htslib/sam.h"

#include <string>
#include <memory>
#include <vector>

namespace gamgee {
/** 
 * @brief simple enum to keep the indices of the genotypes in the PL field of diploid individuals
 */
enum class DiploidPLGenotype { HOM_REF = 0, HET = 1, HOM_VAR = 2};

/**
 * @brief Utility class to manipulate a Variant record.
 */
class Variant {
 public:
  Variant() = default;                                                                                        ///< initializes a null Variant @note this is only used internally by the iterators @warning if you need to create a Variant from scratch, use the builder instead
  explicit Variant(const std::shared_ptr<bcf_hdr_t>& header, const std::shared_ptr<bcf1_t>& body) noexcept;   ///< creates a Variant given htslib objects. @note used by all iterators
  Variant(const Variant& other);                                                                              ///< makes a deep copy of a Variant and it's header. Shared pointers maintain state to all other associated objects correctly.
  Variant(Variant&& other) noexcept;                                                                          ///< moves Variant and it's header accordingly. Shared pointers maintain state to all other associated objects correctly.
  Variant& operator=(const Variant& other);                                                                   ///< deep copy assignment of a Variant and it's header. Shared pointers maintain state to all other associated objects correctly.
  Variant& operator=(Variant&& other) noexcept;                                                               ///< move assignment of a Variant and it's header. Shared pointers maintain state to all other associated objects correctly.

  VariantHeader header() const { return VariantHeader{m_header}; }

  uint32_t chromosome()      const {return uint32_t(m_body->rid);}                                            ///< returns the integer representation of the chromosome. Notice that chromosomes are listed in index order with regards to the header (so a 0-based number). Similar to Picards getReferenceIndex()
  uint32_t alignment_start() const {return uint32_t(m_body->pos+1);}                                          ///< returns a 1-based alignment start position (as you would see in a VCF file). @note the internal encoding is 0-based to mimic that of the BCF files.
  float    qual()            const {return m_body->qual;}                                                     ///< returns the Phred scaled site qual (probability that the site is not reference). See VCF spec.
  uint32_t n_samples()       const {return uint32_t(m_body->n_sample);}                                       ///< returns the number of samples in this Variant record
  uint32_t n_alleles()       const {return uint32_t(m_body->n_allele);}                                       ///< returns the number of alleles in this Variant record

  std::string id() const;                                                                                     ///< returns the variant id field (typically dbsnp id)
  std::string ref() const;                                                                                    ///< returns the ref allele in this Variant record
  std::vector<std::string> alt() const;                                                                       ///< returns the vectors of alt alleles in this Variant record

  // filter field getter
  VariantFilters filters() const;                                                                             ///< returns a vector-like object with all the filters for this record
  bool has_filter(const std::string& filter) const;                                                           ///< checks for the existence of a filter in this record

  // standard individual field getters (a.k.a "format fields")
  IndividualField<IndividualFieldValue<int32_t>> genotype_quals() const;                                            ///< returns a random access object with all the GQ values for all samples contiguously in memory. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<int32_t>> phred_likelihoods() const;                                         ///< returns a random access object with all the PL values for all samples contiguous in memory. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<Genotype> genotypes() const;                                                                   ///< returns a random access object with all the values in a given GT tag for all samples contiguous in memory. @warning Only int8_t GT fields have been tested. @warning Missing GT fields are untested. @warning creates a new object but makes no copies of the underlying values.

  // generic individual field getters (a.k.a "format fields")
  IndividualField<IndividualFieldValue<int32_t>> integer_individual_field(const std::string& tag) const;            ///< returns a random access object with all the values in a give individual field tag in integer format for all samples contiguous in memory.  @warning Only int8_t GT fields have been tested. @warning Missing GT fields are untested. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<float>> float_individual_field(const std::string& tag) const;                ///< returns a random access object with all the values in a give individual field tag in float format for all samples contiguous in memory.  @warning Only int8_t GT fields have been tested. @warning Missing GT fields are untested. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<std::string>> string_individual_field(const std::string& tag) const;         ///< returns a random access object with all the values in a give individual field tag in string format for all samples contiguous in memory. @warning Only int8_t GT fields have been tested. @warning Missing GT fields are untested. @warning creates a new object but makes no copies of the underlying values.

  // generic shared field getters (a.k.a "info fields")
  bool boolean_shared_field(const std::string& tag) const;                    ///< whether or not the tag is present @note bools are treated specially as vector<bool> is impossible given the spec
  SharedField<int32_t> integer_shared_field(const std::string& tag) const;    ///< returns a random access object with all the values in a given shared field tag in integer format contiguous in memory. @warning creates a new vector and copies all underlying values (converting if necessary) to the new vector.
  SharedField<float> float_shared_field(const std::string& tag) const;        ///< returns a random access object with all the values in a given shared field tag in float format for all samples contiguous in memory. @warning creates a new vector and copies all underlying values (converting if necessary) to the new vector.
  SharedField<std::string> string_shared_field(const std::string& tag) const; ///< returns a random access object with all the values in a given shared field tag in string format for all samples contiguous in memory. @warning creates a new vector and copies all underlying values (converting if necessary) to the new vector.

  /**
   * @brief returns a bitset indicating the samples for which the unary predicate is true
   * @tparam any type that can be instantiated in a IndividualFieldIterator<FF>
   */
  template <class FF>                                                                                         
  static boost::dynamic_bitset<> select_if(                                                                   
      const IndividualFieldIterator<FF>& first,         ///< Iterator to the initial position in a sequence. The range includes the element pointed by first.
      const IndividualFieldIterator<FF>& last,          ///< Iterator to the last position in a sequence. The range does not include the element pointed by last.
      const std::function<bool (const FF& value)> pred) ///< Unary predicate function that accepts an element in range [first, last) as argument and returns a value convertible to bool. The value returned indicates whether the element is considered a match in the context of this function. @note The function shall not modify its argument. @note This can either be a function pointer or a function object.
  {
    const auto n_samples = last - first; // HINT: This isn't basic pointer subtraction, it is a customized C++ operator overload that computes the number of samples between two IndividualFieldIterator objects. Took a while to debug when it was not fully tested. -kshakir
    auto selected_samples = boost::dynamic_bitset<>(n_samples);
    auto it = first;
    for (auto i = 0; i < n_samples; i++) {
      selected_samples[i] = pred(*it++);
    }
    return selected_samples;
  }

 private:
  std::shared_ptr<bcf_hdr_t> m_header;                                                                        ///< htslib variant header pointer
  std::shared_ptr<bcf1_t> m_body;                                                                             ///< htslib variant body pointer

  inline bcf_fmt_t* find_individual_field_by_tag(const std::string& tag) const;
  inline bcf_info_t* find_shared_field_by_tag(const std::string& tag) const;
  template <typename TYPE> inline std::vector<TYPE> shared_field(const std::string& tag, const int type) const;

  friend class VariantWriter;
};

}  // end of namespace

#endif /* defined(gamgee__variant__guard) */
