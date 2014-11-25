#ifndef gamgee__variant__guard
#define gamgee__variant__guard

#include "variant_header.h"
#include "individual_field.h"
#include "individual_field_value.h"
#include "shared_field.h"
#include "variant_filters.h"
#include "genotype.h"
#include "utils/variant_utils.h"

#include "htslib/sam.h"
#include "boost/dynamic_bitset.hpp"

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
  Variant& operator=(const Variant& other);                                                                   ///< deep copy assignment of a Variant and it's header. Shared pointers maintain state to all other associated objects correctly.
  Variant(Variant&& other) = default;                                                                         ///< moves Variant and it's header accordingly. Shared pointers maintain state to all other associated objects correctly.
  Variant& operator=(Variant&& other) = default;                                                              ///< move assignment of a Variant and it's header. Shared pointers maintain state to all other associated objects correctly.

  /**
   * @brief returns the header for this variant
   *
   * @note does not deep copy the header; returned VariantHeader object shares existing memory
   */
  VariantHeader header() const {
    // Avoid a deep copy here by constructing using VariantHeader's internal shared header pointer
    return VariantHeader{m_header.m_header};
  }

  bool missing() const { return m_body == nullptr; }                 ///< returns true if this is a default-constructed Variant object with no data

  uint32_t chromosome()         const {return uint32_t(m_body->rid);}                                         ///< returns the integer representation of the chromosome. Notice that chromosomes are listed in index order with regards to the header (so a 0-based number). Similar to Picards getReferenceIndex()
  std::string chromosome_name() const {return header().chromosomes()[chromosome()];}                          ///< returns the name of the chromosome by querying the header.
  uint32_t alignment_start()    const {return uint32_t(m_body->pos+1);}                                       ///< returns a 1-based alignment start position (as you would see in a VCF file). @note the internal encoding is 0-based to mimic that of the BCF files.
  uint32_t alignment_stop()     const {return uint32_t(m_body->pos + m_body->rlen);}                          ///< returns a 1-based alignment stop position, as you would see in a VCF INFO END tag, or the end position of the reference allele if there is no END tag.
  float    qual()               const {return m_body->qual;}                                                  ///< returns the Phred scaled site qual (probability that the site is not reference). See VCF spec.
  uint32_t n_samples()          const {return uint32_t(m_body->n_sample);}                                    ///< returns the number of samples in this Variant record
  uint32_t n_alleles()          const {return uint32_t(m_body->n_allele);}                                    ///< returns the number of alleles in this Variant record including the reference allele

  std::string id() const;                                                                                     ///< returns the variant id field (typically dbsnp id)
  std::string ref() const;                                                                                    ///< returns the ref allele in this Variant record
  std::vector<std::string> alt() const;                                                                       ///< returns the vectors of alt alleles in this Variant record

  // filter field getter
  VariantFilters filters() const;                                                                             ///< returns a vector-like object with all the filters for this record
  bool has_filter(const std::string& filter) const;                                                           ///< checks for the existence of a filter in this record

  // individual field getters (a.k.a "format fields")
  IndividualField<Genotype> genotypes() const;                                                                 ///< special getter for the Genotype (GT) field. Returns a random access object with all the values in a given GT tag for all samples contiguous in memory. @warning Only int8_t GT fields have been tested. @warning Missing GT fields are untested. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<int32_t>> integer_individual_field(const std::string& tag) const;       ///< returns a random access object with all the values in a given individual field tag in integer format for all samples contiguous in memory.  @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<float>> float_individual_field(const std::string& tag) const;           ///< returns a random access object with all the values in a given individual field tag in float format for all samples contiguous in memory.  @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<std::string>> string_individual_field(const std::string& tag) const;    ///< returns a random access object with all the values in a given individual field tag in string format for all samples contiguous in memory. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<int32_t>> individual_field_as_integer(const std::string& tag) const;    ///< same as integer_individual_field but will attempt to convert underlying data to integer if possible. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<float>> individual_field_as_float(const std::string& tag) const;        ///< same as float_individual_field but will attempt to convert underlying data to float if possible. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<std::string>> individual_field_as_string(const std::string& tag) const; ///< same as string_individual_field but will attempt to convert underlying data to string if possible. @warning Only int8_t GT fields have been tested.
  IndividualField<IndividualFieldValue<int32_t>> integer_individual_field(const int32_t index) const;          ///< returns a random access object with all the values in a given individual field tag index in integer format for all samples contiguous in memory. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<float>> float_individual_field(const int32_t index) const;              ///< returns a random access object with all the values in a given individual field tag index in float format for all samples contiguous in memory.   @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<std::string>> string_individual_field(const int32_t index) const;       ///< returns a random access object with all the values in a given individual field tag index in string format for all samples contiguous in memory.   @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<int32_t>> individual_field_as_integer(const int32_t index) const;       ///< same as integer_individual_field but will attempt to convert underlying data to integer if possible. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<float>> individual_field_as_float(const int32_t index) const;           ///< same as float_individual_field but will attempt to convert underlying data to float if possible. @warning creates a new object but makes no copies of the underlying values.
  IndividualField<IndividualFieldValue<std::string>> individual_field_as_string(const int32_t index) const;    ///< same as string_individual_field but will attempt to convert underlying data to string if possible. @warning creates a new object but makes no copies of the underlying values.

  // shared field getters (a.k.a "info fields")
  bool boolean_shared_field(const std::string& tag) const;                       ///< whether or not the tag is present @note bools are treated specially as vector<bool> is impossible given the spec
  SharedField<int32_t> integer_shared_field(const std::string& tag) const;       ///< returns a random access object with all the values in a given shared field tag in integer format contiguous in memory. @warning creates a new object but makes no copies of the underlying values.
  SharedField<float> float_shared_field(const std::string& tag) const;           ///< returns a random access object with all the values in a given shared field tag in float format for all samples contiguous in memory. @warning creates a new object but makes no copies of the underlying values.
  SharedField<std::string> string_shared_field(const std::string& tag) const;    ///< returns a random access object with all the values in a given shared field tag in string format for all samples contiguous in memory. @warning creates a new object but makes no copies of the underlying values.
  SharedField<int32_t> shared_field_as_integer(const std::string& tag) const;    ///< same as integer_shared_field but will attempt to convert underlying data to integer if possible. @warning creates a new object but makes no copies of the underlying values.
  SharedField<float> shared_field_as_float(const std::string& tag) const;        ///< same as float_shared_field but will attempt to convert underlying data to float if possible. @warning creates a new object but makes no copies of the underlying values.
  SharedField<std::string> shared_field_as_string(const std::string& tag) const; ///< same as string_shared_field but will attempt to convert underlying data to string if possible. @warning creates a new object but makes no copies of the underlying values.
  bool boolean_shared_field(const int32_t index) const;                          ///< whether or not the tag with this index is present @note bools are treated specially as vector<bool> is impossible given the spec
  SharedField<int32_t> integer_shared_field(const int32_t index) const;          ///< same as integer_shared_field but will attempt to convert underlying data to integer if possible. @warning creates a new object but makes no copies of the underlying values.
  SharedField<float> float_shared_field(const int32_t index) const;              ///< same as float_shared_field but will attempt to convert underlying data to float if possible. @warning creates a new object but makes no copies of the underlying values.
  SharedField<std::string> string_shared_field(const int32_t index) const;       ///< same as string_shared_field but will attempt to convert underlying data to string if possible. @warning creates a new object but makes no copies of the underlying values.
  SharedField<int32_t> shared_field_as_integer(const int32_t index) const;       ///< same as integer_shared_field but will attempt to convert underlying data to integer if possible. @warning creates a new object but makes no copies of the underlying values.
  SharedField<float> shared_field_as_float(const int32_t index) const;           ///< same as float_shared_field but will attempt to convert underlying data to float if possible. @warning creates a new object but makes no copies of the underlying values.
  SharedField<std::string> shared_field_as_string(const int32_t index) const;    ///< same as string_shared_field but will attempt to convert underlying data to string if possible. @warning creates a new object but makes no copies of the underlying values.

  /**
   * @brief functional-style set logic operations for variant field vectors
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
   * const auto gqs = record.integer_individual_field("GQ"); // a "vector-like" with all the GQs of all samples in this record
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

  /**
   * @brief computes the allele types for all allels (including the reference allele)
   *
   * This function gives you an index vector (AlleleMask) that can be used to
   * query genotypes for snp(), indel(), complex(), ...
   *
   * Complexity is O(N) on the number of alleles. 
   *
   * @return a vector of AlleleType that can be used with many of Genotype member functions
   */
  AlleleMask allele_mask() const;

 private:
  VariantHeader m_header;                                                                        ///< variant header
  std::shared_ptr<bcf1_t> m_body;                                                                ///< htslib variant body pointer

  bcf_fmt_t*  find_individual_field(const std::string& tag) const { return bcf_get_fmt(m_header.m_header.get(), m_body.get(), tag.c_str());  }
  bcf_info_t* find_shared_field(const std::string& tag)     const { return bcf_get_info(m_header.m_header.get(), m_body.get(), tag.c_str()); }
  bcf_fmt_t*  find_individual_field(const uint32_t index) const { return bcf_get_fmt_id(m_body.get(), index); }
  bcf_info_t* find_shared_field(const uint32_t index)     const { return bcf_get_info_id(m_body.get(), index); }
  bool check_field(const int32_t type_field, const int32_t type_value, const int32_t index) const;
  inline AlleleType allele_type_from_difference(const int diff) const;

  template<class FIELD_TYPE, class INDEX_OR_TAG> SharedField<FIELD_TYPE> shared_field_as(const INDEX_OR_TAG& p) const;
  template<class FIELD_TYPE, class INDEX_OR_TAG> IndividualField<IndividualFieldValue<FIELD_TYPE>> individual_field_as(const INDEX_OR_TAG& p) const;

  friend class VariantWriter;
  friend class VariantBuilder; ///< builder needs access to the internals in order to build efficiently

  // TODO: remove this friendship and these mutators after Issue #320 is resolved

  friend class ReferenceBlockSplittingVariantIterator;

  inline void set_alignment_start(const int32_t start) { m_body->pos = start - 1; }
  inline void set_alignment_stop(const int32_t end) { m_body->rlen = end - m_body->pos; }

  inline void set_reference_allele(const char* ref, const int32_t ref_length)
  {
    bcf_unpack(m_body.get(), BCF_UN_STR);
    //Try to avoid calling bcf_update_alleles as it causes significant overhead in
    //terms of memory re-allocation etc.
    //If enough space is available, over-write the current reference allele value
    //Seems like an un-necessary hack, but trust me every single such overhead counts
    if(m_body->rlen >= ref_length)
    {
      memcpy(m_body->d.allele[0], ref, ref_length);	//memcpy is slightly faster than strcpy
      m_body->d.allele[0][ref_length] = '\0';		//append null - why not include in memcpy? see set_reference_allele with char option
      m_body->d.shared_dirty |= BCF1_DIRTY_ALS;
    }
    else
    {
      //expensive - causes significant reallocs within htslib
      //hacks to interpret const ptrs to non-const. Again, saves some overhead
      //No need to free d.allele[0] because it points to shared string within bcf1_1.shared
      m_body->d.allele[0] = const_cast<char*>(ref);
      //Re-use same d.allele char** as update_alleles writes to different region of memory anyway
      bcf_update_alleles(const_cast<const bcf_hdr_t*>(m_header.m_header.get()), m_body.get(), const_cast<const char**>(m_body->d.allele), m_body->n_allele);
    }
  }
  inline void set_reference_allele(const char* ref)  { set_reference_allele(ref, static_cast<int32_t>(strlen(ref))); }
  inline void set_reference_allele(const char ref_base)  { set_reference_allele(&ref_base, 1); }
};

}  // end of namespace

#endif /* defined(gamgee__variant__guard) */
