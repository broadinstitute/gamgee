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
  uint32_t alignment_stop()  const {return uint32_t(m_body->pos + m_body->rlen);}                             ///< returns a 1-based alignment stop position, as you would see in a VCF INFO END tag, or the end position of the reference allele if there is no END tag.
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
  IndividualField<IndividualFieldValue<int32_t>> genotype_quals() const;    ///< returns a random access object with all the GQ values for all samples contiguously in memory. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<int32_t>> phred_likelihoods() const; ///< returns a random access object with all the PL values for all samples contiguous in memory. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<Genotype> genotypes() const;                              ///< returns a random access object with all the values in a given GT tag for all samples contiguous in memory. @warning Only int8_t GT fields have been tested. @warning Missing GT fields are untested. @warning creates a new object but makes no copies of the underlying values.

  // generic individual field getters (a.k.a "format fields")
  IndividualField<IndividualFieldValue<int32_t>> integer_individual_field(const std::string& tag) const;    ///< returns a random access object with all the values in a given individual field tag in integer format for all samples contiguous in memory.  @warning Only int8_t GT fields have been tested. @warning Missing GT fields are untested. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<float>> float_individual_field(const std::string& tag) const;        ///< returns a random access object with all the values in a given individual field tag in float format for all samples contiguous in memory.  @warning Only int8_t GT fields have been tested. @warning Missing GT fields are untested. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<std::string>> string_individual_field(const std::string& tag) const; ///< returns a random access object with all the values in a given individual field tag in string format for all samples contiguous in memory. @warning Only int8_t GT fields have been tested. @warning Missing GT fields are untested. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<int32_t>> individual_field_as_integer(const std::string& tag) const;    ///< same as integer_format_field but will attempt to convert underlying data to integer if possible. @warning Only int8_t GT fields have been tested. @warning Missing GT fields are untested. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<float>> individual_field_as_float(const std::string& tag) const;        ///< same as float_format_field but will attempt to convert underlying data to float if possible. @warning Only int8_t GT fields have been tested. @warning Missing GT fields are untested. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<std::string>> individual_field_as_string(const std::string& tag) const; ///< same as string_format_field but will attempt to convert underlying data to string if possible. @warning Only int8_t GT fields have been tested. @warning Missing GT fields are untested. @warning creates a new object but makes no copies of the underlying values.


  // generic shared field getters (a.k.a "info fields")
  bool boolean_shared_field(const std::string& tag) const;                    ///< whether or not the tag is present @note bools are treated specially as vector<bool> is impossible given the spec
  SharedField<int32_t> integer_shared_field(const std::string& tag) const;    ///< returns a random access object with all the values in a given shared field tag in integer format contiguous in memory. @warning creates a new object but makes no copies of the underlying values.
  SharedField<float> float_shared_field(const std::string& tag) const;        ///< returns a random access object with all the values in a given shared field tag in float format for all samples contiguous in memory. @warning creates a new object but makes no copies of the underlying values.
  SharedField<std::string> string_shared_field(const std::string& tag) const; ///< returns a random access object with all the values in a given shared field tag in string format for all samples contiguous in memory. @warning creates a new object but makes no copies of the underlying values.
  SharedField<int32_t> shared_field_as_integer(const std::string& tag) const;    ///< same as integer_shared_field but will attempt to convert underlying data to integer if possible. @warning creates a new object but makes no copies of the underlying values.
  SharedField<float> shared_field_as_float(const std::string& tag) const;        ///< same as float_shared_field but will attempt to convert underlying data to float if possible. @warning creates a new object but makes no copies of the underlying values.
  SharedField<std::string> shared_field_as_string(const std::string& tag) const; ///< same as string_shared_field but will attempt to convert underlying data to string if possible. @warning creates a new object but makes no copies of the underlying values.
  /**
   * @brief functional-stlye set logic operations for variant field vectors
   *
   * This function applies the unary predicate pred to every element between 
   * first and last in the container without modifying them. It then produces 
   * a bitset with as many elements as the full iteration (from first to last)
   * For all elements where the predicate returns true, the corresponding bit 
   * in the bitset will be 1 (set). For all elements where the predicate returns 
   * false, the corresponding bit will be 0 (unset)
   * 
   * The great advantage of this functionality is that you can perform one 
   * operation across all elements of an iterator benefiting from data locality
   * and cache prefetching which gets translated into a more manageable format
   * (bitset) for subsequent set-logic transformations. This will be particularly
   * fast for very large datasets, for example, in a VCF file with 1000+ samples,
   * running several select_if's end to end on all the samples and then combining
   * the resulting bitsets using set-logic, is much faster than having one 
   * single iteration updating several IndividualFields due to it's distance in 
   * the computer's memory.
   *
   *
   * For example:
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * const auto genotypes = record.genotypes(); // a "vector-like" with the genotypes of all samples in this record
   * const auto gqs = record.genotype_quals(); // a "vector-like" with all the GQs of all samples in this record
   * const auto hets = Variant::select_if(genotypes.begin(), genotypes.end(), [](const auto& g) { return g.het(); }); // returns a bit set with all hets marked with 1's
   * const auto pass_gqs = Variant::select_if(gqs.begin(), gqs.end(), [](const auto& gq) { return gq[0] > 20; }); // returns a bit set with every sample with gq > 20 marked with 1's
   * const auto high_qual_hets = hets & pass_gqs; // a bit set with all the samples that are het and have gq > 20
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * in the above example you see two invocations of select_if using the same
   * ITER but different VALUEs but the user is oblivious to that as he simply
   * defines the type in his lambda function as const auto&
   * 
   * This function is declared with template template arguments to facilitate
   * parameter type deduction and enable usage without specifying the actual
   * types being used. The whole goal here is to keep the user happily using
   * auto everywhere and unaware of the iterator and value types underlying the
   * data structures.
   *
   * @note pred can either be a function pointer, a function object or a lambda function. 
   * @note This function can be called directly (ignoring the template parameters) as all the template parameters can be deduced from the function parameters. 
   * @tparam ITER any iterator that has operator- defined to return the difference in number of elements between two ITER iterators
   * @tparam VALUE the class of the objects ITER is iterating over. (e.g. in IndividualFieldIterator<Genotype> Genotype is the VALUE, IndividualFieldIterator is the ITER
   * @param first iterator to the initial position in a sequence. The range includes the element pointed by first.
   * @param last iterator to the last position in a sequence. The range does not include the element pointed by last.
   * @param pred unary predicate (lambda) function that accepts an element in range [first, last) as argument and returns a value convertible to bool. The value returned indicates whether the element is considered a match in the context of this function. 
   * @return a bitset indicating the samples for which the unary predicate is true
   */
  template <class VALUE, template<class> class ITER>                                                                                         
  static boost::dynamic_bitset<> select_if(                                                                   
      const ITER<VALUE>& first,         
      const ITER<VALUE>& last,          
      const std::function<bool (const decltype(*first)& value)> pred)  // decltype(*first) will always yield an l-value reference to VALUE (VALUE&). I just added the & here to be explicit about my intention of capturing the parameter by reference.
  {
    const auto n_samples = last - first; 
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
