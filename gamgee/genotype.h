#ifndef gamgee__genotype__guard
#define gamgee__genotype__guard

#include "utils/variant_utils.h"

#include "utils/hts_memory.h"
#include "utils/utils.h"
#include "utils/variant_field_type.h"
#include "utils/genotype_utils.h"

#include <memory>
#include <utility>
#include <stdexcept>

namespace gamgee {

/**
 * @brief Encodes a genotype.
 */
class Genotype{

 public:
  /**
   * @brief Constructs a genotype.
   * @param body The shared memory variant "line" from a vcf, or bcf.
   * @param format_ptr The GT field from the line.
   * @param data_ptr The GT for this sample.
   */
  Genotype(const std::shared_ptr<bcf1_t>& body, const bcf_fmt_t* const format_ptr, const uint8_t* data_ptr);

  /**
   * @brief copying of the Genotype object is not allowed.
   */
  Genotype(const Genotype& other) = delete;

  /**
   * @brief copying of the Genotype object is not allowed.
   * @param other Other genotype.
   */
  Genotype& operator=(const Genotype& other) = delete;

  /**
   * @brief Explicit default as recommended by many threads on stackoverflow.
   * @param other Other genotype.
   */
  Genotype(Genotype&& other) = default;

  /**
   * @brief Explicit default as recommended by many threads on stackoverflow.
   * @param other Other genotype.
   */
  Genotype& operator=(Genotype&& other) = default;

  /**
   * @brief Explicit default as recommended by many threads on stackoverflow.
   */
  ~Genotype() = default;

  /**
   * @brief Checks if another genotype does not equal this genotype.
   * @param other The other genotype to compare.
   * @return False if other genotype equals this genotype.
   */
  bool operator!=(const Genotype& other) const;

  /**
   * @brief Checks if another genotype equals this genotype.
   * @param other The other genotype to compare.
   * @return True if other genotype equals this genotype.
   * @note No string comparison is done. This operation is supposed to be fast.
   */
  bool operator==(const Genotype& other) const;

  // only for diploids
  /**
   * @brief Checks if this genotype vector is any type of heterozygous call.
   * @return True if this GT is a het.
   * @note only for diploids.
   */
  bool het() const;

  /**
   * @brief Checks if this genotype vector is a heterozygous call and none of the alleles is the reference.
   * @return True if this GT is a het and none of the alleles is the reference.
   * @note only for diploids.
   */
  bool non_ref_het() const;

  /**
   * @brief A bit encoding for the first two alleles.
   * @return A bit encoding for the first two alleles.
   * @note only for diploids, returns false otherwise.
   */
  uint32_t fast_diploid_key_generation() const;

  //for all ploidies
  /**
   * @brief Checks if this genotype vector is a homozygous call that is non-reference.
   * @return True if this GT is a hom var.
   */
  bool hom_var() const;

  /**
   * @brief Checks if this genotype vector is a homozygous call that is reference.
   * @return True if this GT is a hom ref.
   */
  bool hom_ref() const;

  /**
   * @brief Checks if all alleles are missing.
   * @return True if all alleles are missing.
   * @warning Missing GT fields are untested.
   */
  bool missing() const;

  // referencing alleles inside a sample's genotype

  /**
   * @brief Returns a vector with all the allele strings
   * @return a vector with all the allele strings
   */
  std::vector<std::string> allele_strings() const;

  /**
   * @brief Returns a vector with all the allele keys
   * @return a vector with all the allele keys
   */
  std::vector<int32_t> allele_keys() const;

  /**
   * @brief Returns the allele string at index
   * @param index Zero based allele index
   * @return the allele string
   */
  std::string allele_string(const uint32_t index) const;

  /**
   * @brief Returns the allele key within this line
   * @param index Zero based allele index
   * @return the allele key within this line
   */
  int32_t allele_key(const uint32_t index) const;

  /**
   * @brief Returns the allele key within this line
   * @param index Zero based allele index
   * @return the allele key within this line
   */
  int32_t operator[](const uint32_t index) const;

  /**
   * @brief Returns the number of alleles
   * @return the number of alleles
   */
  uint32_t size() const;

  /**
   * @brief whether or not this genotype represents a snp
   *
   * In the true sense of the word, single nucleotide polymorphism restricts
   * the number of loci (nucleotides) that are polymorphic, but not the number
   * of potentially different alleles it may have (as long as they are all
   * single nucleotide polymorphisms). This function will check that at least
   * one of these SNPs exists and that no other type of allele (insertion or
   * deletion) is present. Multiple SNPs will return true.
   *
   * @warning complexity is O(n) in the number of alleles. If you are using
   * many of these convenience functions (snp(), insertion(), deletion(),
   * indel(), complex(),...), you will be better off implementing one loop that
   * makes all the checks in one pass instead of calling many O(n) functions.
   *
   * @return whether or not there is at least one snp in this genotype and
   * nothing else but snps and reference alleles
   *
   */
  bool snp(const AlleleMask& mask) const;

  /** 
   * @brief whether or not this genotype represents an insertion
   *
   * This function will check that at least one insertion exists and that
   * no other type of allele (snp or deletion) is present. Multiple insertions
   * will return true.
   *
   * @warning complexity is O(n) in the number of alleles. If you are using
   * many of these convenience functions (snp(), insertion(), deletion(),
   * indel(), complex(),...), you will be better off implementing one loop that
   * makes all the checks in one pass instead of calling many O(n) functions.
   *
   * @return whether or not there is at least one insertion in this genotype and
   * nothing else but insertions and reference alleles
   */
  bool insertion(const AlleleMask& mask) const;

  /** 
   * @brief whether or not this genotype represents an deletion
   *
   * This function will check that at least one deletion exists and that
   * no other type of allele (snp or insertion) is present. Multiple deletions
   * will return true.
   *
   * @warning complexity is O(n) in the number of alleles. If you are using
   * many of these convenience functions (snp(), insertion(), deletion(),
   * indel(), complex(),...), you will be better off implementing one loop that
   * makes all the checks in one pass instead of calling many O(n) functions.
   *
   * @return whether or not there is at least one deletion in this genotype and
   * nothing else but deletions and reference alleles
   */
  bool deletion(const AlleleMask& mask) const;

  /** 
   * @brief whether or not this genotype represents an insertion or deletion
   *
   * This function will check that at least one insertion or deletion exists
   * and that no other type of allele (snp or insertion) is present. Multiple
   * deletions will return true.
   *
   * @note this is not the same as insertion(mask) || deletion(mask) because it
   * will also tolerate sites with insertions and deletions, while both other
   * functions would return false to such a site.
   *
   * @warning complexity is O(n) in the number of alleles. If you are using
   * many of these convenience functions (snp(), insertion(), deletion(),
   * indel(), complex(),...), you will be better off implementing one loop that
   * makes all the checks in one pass instead of calling many O(n) functions.
   *
   * @return whether or not there is at least one insertion or deletion in this
   * genotype and nothing else but insertions, deletions and reference alleles
   */
  bool indel(const AlleleMask& mask) const;

  /** 
   * @brief whether or not this genotype has *at most* one alternate allele
   *
   * This function will check whether this is a simple heterozygous or 
   * homozygous site where both alleles are either the same, or reference and 
   * one alternate allele. No two different alleles would pass.
   *
   * @warning complexity is O(n) in the number of alleles. If you are using
   * many of these convenience functions (snp(), insertion(), deletion(),
   * indel(), complex(),...), you will be better off implementing one loop that
   * makes all the checks in one pass instead of calling many O(n) functions.
   *
   * @return whether or not all alleles are either the same, or there is one alt allele 
   * mixed with reference alleles.
   */
  bool biallelic() const;

  /**
   * @brief literally the negation of biallelic(mask)
   *
   * @warning complexity is O(n) in the number of alleles. If you are using
   * many of these convenience functions (snp(), insertion(), deletion(),
   * indel(), complex(),...), you will be better off implementing one loop that
   * makes all the checks in one pass instead of calling many O(n) functions.
   *
   * @return whether or not there are more than one alt allele in this record
   */
  bool complex() const { return !biallelic(); }

  /**
   * @brief identifies variants with two different types of alleles 
   *
   * @warning complexity is O(n) in the number of alleles. If you are using
   * many of these convenience functions (snp(), insertion(), deletion(),
   * indel(), complex(),...), you will be better off implementing one loop that
   * makes all the checks in one pass instead of calling many O(n) functions.
   *
   * @return whether or not there are more than one types of alt allele in this record
   */
  bool mixed() const;

  bool variant() const {
    return !missing() && !hom_ref();
  }

  /**
   * @brief Converts a vector of allele indices representing a genotype into BCF-encoded
   *        format suitable for passing to VariantBuilder::set_genotypes(). No phasing
   *        is added.
   *
   *        Example: if you want to encode the genotype 0/1, create a vector with {0, 1}
   *                 and then pass it to this function
   */
  static inline void encode_genotype(std::vector<int32_t>& alleles) {
    encode_genotype(alleles, false);
  }

  /**
   * @brief Converts a vector of allele indices representing a genotype into BCF-encoded
   *        format suitable for passing to VariantBuilder::set_genotypes(), and also
   *        allows you to phase all alleles
   *
   *        Example: if you want to encode the genotype 0|1, create a vector with {0, 1}
   *                 and then pass it to this function with phase_all_alleles set to true
   */
  static inline void encode_genotype(std::vector<int32_t>& alleles, bool phase_all_alleles) {
    for ( auto allele_index = 0u; allele_index < alleles.size(); ++allele_index ) {
      // Only legal value below -1 is the int32 vector end value
      if ( alleles[allele_index] < -1 && alleles[allele_index] != bcf_int32_vector_end ) {
        throw std::invalid_argument{"Genotype vector must consist only of allele indices, -1 for missing values, or vector end values"};
      }
      // Do not modify vector end values
      else if ( alleles[allele_index] != bcf_int32_vector_end ) {
        alleles[allele_index] = (alleles[allele_index] + 1) << 1 | (phase_all_alleles && allele_index > 0u ? 1 : 0);
      }
    }
  }

  /**
   * @brief Converts multiple vectors of allele indices representing genotypes into
   *        BCF_encoded format suitable for passing to VariantBuilder::set_genotypes().
   *        No phasing is added.
   *
   *        Example: if you want to encode the genotypes 0/1 and 1/1, create a vector
   *                 with { {0, 1}, {1, 1} } and pass it to this function
   */
  static inline void encode_genotypes(std::vector<std::vector<int32_t>>& multiple_genotypes) {
    for ( auto& genotype : multiple_genotypes ) {
      encode_genotype(genotype, false);
    }
  }

 private:
  std::shared_ptr<bcf1_t> m_body;
  const bcf_fmt_t* m_format_ptr;
  const uint8_t* m_data_ptr;

  bool allele_is_type_or_ref(const AlleleType& type, const std::vector<int32_t>& keys, const AlleleMask& mask) const;
};

}

#endif
