#ifndef gamgee__variant__guard
#define gamgee__variant__guard

#include "variant_header.h"
#include "variant_field.h"
#include "variant_field_value.h"
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
  Variant() = default;                                                                                                ///< @brief initializes a null Variant @note this is only used internally by the iterators @warning if you need to create a Variant from scratch, use the builder instead
  explicit Variant(const std::shared_ptr<bcf_hdr_t>& header, const std::shared_ptr<bcf1_t>& body) noexcept;           ///< @brief creates a Variant given htslib objects. @note used by all iterators
  Variant(const Variant& other);                                                                                      ///< @brief makes a deep copy of a Variant and it's header. Shared pointers maintain state to all other associated objects correctly.
  Variant(Variant&& other) noexcept;                                                                                  ///< @brief moves Variant and it's header accordingly. Shared pointers maintain state to all other associated objects correctly.
  Variant& operator=(const Variant& other);                                                                           ///< @brief deep copy assignment of a Variant and it's header. Shared pointers maintain state to all other associated objects correctly.
  Variant& operator=(Variant&& other) noexcept;                                                                       ///< @brief move assignment of a Variant and it's header. Shared pointers maintain state to all other associated objects correctly.

  VariantHeader header() const { return VariantHeader{m_header}; }

  uint32_t chromosome()      const {return uint32_t(m_body->rid);}                                                    ///< @brief returns the integer representation of the chromosome. Notice that chromosomes are listed in index order with regards to the header (so a 0-based number). Similar to Picards getReferenceIndex()
  uint32_t alignment_start() const {return uint32_t(m_body->pos+1);}                                                  ///< @brief returns a 1-based alignment start position (as you would see in a VCF file). @note the internal encoding is 0-based to mimic that of the BCF files.
  float    qual()            const {return m_body->qual;}                                                             ///< @brief returns the Phred scaled site qual (probability that the site is not reference). See VCF spec.
  uint32_t n_samples()       const {return uint32_t(m_body->n_sample);}                                               ///< @brief returns the number of samples in this Variant record
  uint32_t n_alleles()       const {return uint32_t(m_body->n_allele);}                                               ///< @brief returns the number of alleles in this Variant record

  std::string id() const;                                                                                             ///< @brief returns the variant id field (typically dbsnp id)
  std::string ref() const;                                                                                            ///< @brief returns the ref allele in this Variant record
  std::vector<std::string> alt() const;                                                                               ///< @brief returns the vectors of alt alleles in this Variant record

  // filter field getter
  VariantFilters filters() const;                                                                                     ///< @brief returns a vector-like object with all the filters for this record
  bool has_filter(const std::string& filter) const;                                                                   ///< @brief checks for the existence of a filter in this record

  // standard format field getters
  bool is_hom_ref(const uint32_t sample_index) const;                                                                 ///< @brief whether or not the sample in sample_index is hom_ref @warning current implementation is solely based on PL (this will be changed in the future)
  bool is_het(const uint32_t sample_index) const;                                                                     ///< @brief whether or not the sample in sample_index is het @warning current implementation is solely based on PL (this will be changed in the future)
  bool is_hom_var(const uint32_t sample_index) const;                                                                 ///< @brief whether or not the sample in sample_index is hom_var @warning current implementation is solely based on PL (this will be changed in the future)
  VariantField<VariantFieldValue<int32_t>> genotype_quals() const;                                                    ///< @brief returns a random access object with all the GQ values for all samples contiguously in memory.
  VariantField<VariantFieldValue<int32_t>> phred_likelihoods() const;                                                 ///< @brief returns a random access object with all the PL values for all samples contiguous in memory.

  // generic format field getters
  VariantField<VariantFieldValue<int32_t>> generic_integer_format_field(const std::string& tag) const;                ///< @brief returns a random access object with all the values in a given format field tag in integer format for all samples contiguous in memory.
  VariantField<VariantFieldValue<float>> generic_float_format_field(const std::string& tag) const;                    ///< @brief returns a random access object with all the values in a given format field tag in float format for all samples contiguous in memory.
  VariantField<VariantFieldValue<std::string>> generic_string_format_field(const std::string& tag) const;             ///< @brief returns a random access object with all the values in a given format field tag in string format for all samples contiguous in memory. @warning not working at the moment due to bug in htslib

  VariantField<Genotype> genotypes() const;                                                                           ///< @brief returns a random access object with all the values in a given GT tag for all samples contiguous in memory. @warning Only int8_t GT fields have been tested. @warning Missing GT fields are untested.

  // generic info field getters
  // NOTE: Errors in the return info field type are currently returned as empty vectors, not exceptions.
  std::vector<int32_t> generic_integer_info_field(const std::string& tag) const;                                      ///< @brief returns a random access object with all the values in a given info field tag in integer format for all samples contiguous in memory.
  std::vector<float> generic_float_info_field(const std::string& tag) const;                                          ///< @brief returns a random access object with all the values in a given info field tag in float format for all samples contiguous in memory.
  std::vector<std::string> generic_string_info_field(const std::string& tag) const;                                   ///< @brief returns a random access object with all the values in a given info field tag in string format for all samples contiguous in memory.
  std::vector<bool> generic_boolean_info_field(const std::string& tag) const;                                         ///< @brief returns a random access object with all the values in a given info field tag in boolean format for all samples contiguous in memory.

  /**
  * @brief returns a bitset indicating the samples selected according to unary predicate.
  * */
  template <class FF>
  static boost::dynamic_bitset<> select_if(
      const VariantFieldIterator<FF>& first,                                                                          ///< Iterator to the initial position in a sequence. The range includes the element pointed by first.
      const VariantFieldIterator<FF>& last,                                                                           ///< Iterator to the last position in a sequence. The range does not include the element pointed by last.
      const std::function<bool (const FF& value)> pred)                                                               ///< Unary predicate function that accepts an element in range [first, last) as argument and returns a value convertible to bool. The value returned indicates whether the element is considered a match in the context of this function. @note The function shall not modify its argument. @note This can either be a function pointer or a function object.
  {
    const auto n_samples = last - first; // HINT: This isn't basic pointer subtraction, it is a customized C++ operator overload that computes the number of samples between two VariantFieldIterator objects. Took a while to debug when it was not fully tested. -kshakir
    auto selected_samples = boost::dynamic_bitset<>(n_samples);
    auto it = first;
    for (auto i = 0; i < n_samples; i++) {
      selected_samples[i] = pred(*it++);
    }
    return selected_samples;
  }

 private:
  std::shared_ptr<bcf_hdr_t> m_header;                                                                                ///< @brief htslib variant header pointer
  std::shared_ptr<bcf1_t> m_body;                                                                                     ///< @brief htslib variant body pointer

  inline bcf_fmt_t* find_format_field_by_tag(const std::string& tag) const;
  template <typename TYPE> inline std::vector<TYPE> generic_info_field(const std::string& tag, const int type) const;
  inline bool is_this_genotype(const DiploidPLGenotype& genotype, const uint32_t sample_index) const;

  friend class VariantWriter;
};

}  // end of namespace

#endif /* defined(gamgee__variant__guard) */
