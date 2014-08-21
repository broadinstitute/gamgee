#ifndef gamgee__genotype__guard
#define gamgee__genotype__guard

#include "utils/hts_memory.h"
#include "utils/utils.h"
#include "utils/variant_field_type.h"
#include "utils/genotype_utils.h"

#include <memory>
#include <utility>

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
   * @brief Checks if this genotype vector is a homozygous call that is non-reference.
   * @return True if this GT is a hom var.
   * @note only for diploids, returns false otherwise.
   */
  bool hom_var() const;
  /**
   * @brief Checks if this genotype vector is a homozygous call that is reference.
   * @return True if this GT is a hom ref.
   * @note only for diploids, returns false otherwise.
   */
  bool hom_ref() const;
  /**
   * @brief A bit encoding for the first two alleles.
   * @return A bit encoding for the first two alleles.
   * @note only for diploids, returns false otherwise.
   */
  uint32_t fast_diploid_key_generation() const;

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
  std::vector<std::string> alleles_strings() const;

  /**
   * @brief Returns a vector with all the allele keys
   * @return a vector with all the allele keys
   */
  std::vector<int32_t> alleles_keys() const;

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

  // TODO: Discussions around having VariantFieldValueIterator work on Genotype objects too, or return some other iterator.
  // ? begin() const;
  // ? end() const;

  private:
  const std::shared_ptr<bcf1_t> m_body;
  const bcf_fmt_t* const m_format_ptr;
  const uint8_t* m_data_ptr;
};

}

#endif
