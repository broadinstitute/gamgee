#ifndef gamgee__variant_builder__guard
#define gamgee__variant_builder__guard

#include "variant.h"
#include "variant_builder_shared_region.h"
#include "variant_builder_individual_region.h"

#include <vector>
#include <string>

namespace gamgee {

/**
 * @brief Internal VariantBuilder class to represent "core" (non-data-region) fields such as the alignment start and the qual.
 *
 * Necessary because we need to distinguish between the case where the user explicitly sets
 * one of these to a missing value vs. the case where the user doesn't touch a field at all.
 */
template<class FIELD_TYPE>
class VariantBuilderCoreField {
 public:
  VariantBuilderCoreField() : m_field_value{}, m_is_set{false} {}

  void set(const FIELD_TYPE value) { m_field_value = value; m_is_set = true; }
  FIELD_TYPE field_value() const { return m_field_value; }
  bool is_set() const { return m_is_set; }
  void clear() { m_is_set = false; }

 private:
  FIELD_TYPE m_field_value;
  bool m_is_set;
};


/**
 * @brief VariantBuilder: construct Variant records from scratch
 *                        (and, coming soon, from existing Variant records)
 *
 * To use, first create a VariantHeader appropriate for the final
 * record you intend to create (with all shared/individual fields,
 * samples, and contigs pre-declared), then use that header to
 * instantiate a builder:
 *
 * auto builder = VariantBuilder{header};
 *
 * You should create ONE builder per file you intend to output, and re-use
 * it across records, calling clear() in between each record. Do NOT
 * create a new builder for each record -- creating/destroying builders is
 * an expensive process involving many memory allocations and deallocations
 * of internal data structures as well as costly header lookups. Ie.,
 * the correct way to use a builder is:
 *
 * auto builder = VariantBuilder{header};
 * for ( each record you want to create ) {
 *   // use existing builder to build new record
 *   builder.clear();
 * }
 *
 * Once you have a builder you can call setter functions in a chained fashion
 * as follows:
 *
 * auto variant = builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A")
 *                       .set_alt_alleles({"C", "T"}).set_genotypes(std::move(my_genotypes))
 *                       .build();
 *
 * See the comments further down in this header for instructions and tips on using the various
 * kinds of setters.
 *
 * You must at a minimum set the required chromosome, alignment start, and ref allele fields
 * (unless you have disabled validation, but then you will just get an invalid Variant record).
 *
 * The VariantBuilder API is designed to allow you to be efficient when you want to be (eg.,
 * moving existing data into the builder, setting individual fields in bulk rather than by
 * sample, etc.), and lazy when you don't care about efficiency. In general, the more efficient
 * API functions require a bit more work to use than the less efficient ones. See the discussions
 * below of the efficiency of the various options available to you for setting fields.
 *
 * Some general guidelines:
 *
 * -For removing fields you generally have the option of either calling the appropriate remove_*
 *  API function, OR passing in a missing/empty value to the appropriate set_* function. Both
 *  options are equivalent and will result in the field being removed.
 *
 * -For individual fields, it is an error to request both bulk and per-sample changes to the
 *  same field (without calling clear() in between). This is because it would be too costly
 *  to reconcile the two kinds of changes.
 *
 * -Setting by field/sample index is more efficient than setting by field/sample names, provided
 *  you look up the indices ONCE and cache them.
 *
 * -Setting by move/r-value is more efficient than setting by l-value, as the functions that take
 *  an l-value reference assume that you want to keep your data and therefore make a copy of it.
 *
 * -Setting individual fields in bulk is more efficient than setting by sample, and within the bulk
 *  setters the functions that take a single-dimensional vector are more efficient (but harder to
 *  use) than the functions that take a two-dimensional vector.
 *
 * -Disabling validation is possible, and will certainly improve performance, however if you
 *  attempt to perform an action that would have been prevented by validation checks (such as
 *  setting a non-existent field) you WILL get undefined behavior in your program. You should
 *  ONLY disable validation if you're extremely confident that your program will not take any
 *  invalid/incorrect actions, and that the data you pass to the builder will always be valid.
 */
class VariantBuilder {
 public:
  explicit VariantBuilder(const VariantHeader& header);

  // VariantBuilders are moveable but not copyable, and use default destruction
  VariantBuilder(VariantBuilder&& other) = default;
  VariantBuilder& operator=(VariantBuilder&& other) = default;
  VariantBuilder(const VariantBuilder& other) = delete;
  VariantBuilder& operator=(const VariantBuilder& other) = delete;
  ~VariantBuilder() = default;

  VariantHeader header() const { return VariantHeader{m_header}; }

  /**
   * @brief Disable or enable all validation checks in this builder.
   *
   * Passing in false disables validation, passing in true enables it.
   * Validation is on by default in new builders.
   *
   * @warning Disabling validation will improve performance at the cost of safety.
   *          Only do so if you're extremely confident that your program will not
   *          take any invalid/incorrect actions, and that the data you pass to the
   *          builder will always be valid.
   */
  VariantBuilder& set_enable_validation(const bool enable_validation);

  /******************************************************************************
   *
   * Functions for setting core site-level fields
   *
   ******************************************************************************/

  VariantBuilder& set_chromosome(const uint32_t chromosome);
  VariantBuilder& set_chromosome(const std::string& chromosome);
  VariantBuilder& set_alignment_start(const uint32_t alignment_start);
  VariantBuilder& set_alignment_stop(const uint32_t alignment_stop);
  VariantBuilder& set_qual(const float qual);
  VariantBuilder& set_id(const std::string& id);
  VariantBuilder& set_ref_allele(const std::string& ref_allele);
  VariantBuilder& set_alt_allele(const std::string& alt_allele);
  VariantBuilder& set_alt_alleles(const std::vector<std::string>& alt_alleles);
  VariantBuilder& set_filters(const VariantFilters& filters);
  VariantBuilder& set_filters(const std::vector<std::string>& filters);
  VariantBuilder& set_filters(const std::vector<int32_t>& filters);

  /******************************************************************************
   *
   * Functions for removing core site-level fields
   * (Note that we don't have removal functions for the required fields)
   *
   ******************************************************************************/

  VariantBuilder& remove_alignment_stop();
  VariantBuilder& remove_qual();
  VariantBuilder& remove_id();
  VariantBuilder& remove_alt_alleles();
  VariantBuilder& remove_filters();

  /******************************************************************************
   *
   * Functions for setting shared/INFO fields
   *
   ******************************************************************************/

  VariantBuilder& set_integer_shared_field(const std::string& tag, const int32_t value);
  VariantBuilder& set_integer_shared_field(const std::string& tag, const std::vector<int32_t>& values);
  VariantBuilder& set_integer_shared_field(const uint32_t index, const int32_t value);
  VariantBuilder& set_integer_shared_field(const uint32_t index, const std::vector<int32_t>& values);
  VariantBuilder& set_float_shared_field(const std::string& tag, const float value);
  VariantBuilder& set_float_shared_field(const std::string& tag, const std::vector<float>& values);
  VariantBuilder& set_float_shared_field(const uint32_t index, const float value);
  VariantBuilder& set_float_shared_field(const uint32_t index, const std::vector<float>& values);
  VariantBuilder& set_string_shared_field(const std::string& tag, const std::string& value);
  VariantBuilder& set_string_shared_field(const uint32_t index, const std::string& value);
  VariantBuilder& set_boolean_shared_field(const std::string& tag);
  VariantBuilder& set_boolean_shared_field(const uint32_t index);

  /******************************************************************************
   *
   * Functions for removing shared/INFO fields
   *
   ******************************************************************************/

  VariantBuilder& remove_shared_field(const std::string& tag);
  VariantBuilder& remove_shared_field(const uint32_t field_index);
  VariantBuilder& remove_shared_fields(const std::vector<std::string>& tags);
  VariantBuilder& remove_shared_fields(const std::vector<uint32_t>& field_indices);

  /******************************************************************************
   *
   * Functions for setting individual/FORMAT fields in bulk
   * (ie., setting all values for all samples at once)
   *
   * -These are much more efficient than the per-sample setters.
   *
   * -Setting by field index is more efficient than setting by field name,
   *  provided that you look up the indices ONCE in the header and cache
   *  them instead of looking them up for every record.
   *
   * -The bulk setting functions that take a vector by move/r-value reference
   *  are much more efficient than the bulk setting functions that take a vector
   *  by l-value reference, since the latter make a copy of the data.
   *
   * -The bulk setting functions that take a flattened (one-dimensional) vector
   *  are much more efficient than the functions that take a two-dimensional vector,
   *  but the flattened approach requires the user to pad sample values with
   *  missing/end of vector values manually, whereas with the two-dimensional
   *  approach the user does not need to worry about padding sample values.
   *
   *  For example, with 4 samples and an integer individual field with a
   *  varied number of values per sample, you could pass in the following
   *  nested vector:
   *
   *  { {1, 2}, {3}, {}, {5, 6, 7} }
   *
   *  with each inner vector representing the values for one sample.
   *
   *  However, to use the flattened vector API functions you'd have to pass in
   *  this equivalent single-dimensional vector:
   *
   *  { 1, 2, int32_vector_end,
   *    3, int32_vector_end, int32_vector_end,
   *    int32_missing_value, int32_vector_end, int32_vector_end,
   *    5, 6, 7 }
   *
   *  The advantage of the flattened approach is much greater efficiency
   *  in terms of data locality and memory allocations.
   *
   * -The Genotype setter functions require you to encode your genotype data
   *  with one of the Genotype::encode_genotype() functions before passing it in.
   *
   *  For example, if you want to encode the genotype 0/1, create a vector with
   *  {0, 1}, pass it to Genotype::encode_genotype(), then pass it to the appropriate
   *  set_genotypes() function below.
   *
   ******************************************************************************/

  VariantBuilder& set_genotypes(const std::vector<int32_t>& genotypes_for_all_samples);
  VariantBuilder& set_genotypes(std::vector<int32_t>&& genotypes_for_all_samples);
  VariantBuilder& set_genotypes(const std::vector<std::vector<int32_t>>& genotypes_for_all_samples);
  VariantBuilder& set_genotypes(std::vector<std::vector<int32_t>>&& genotypes_for_all_samples);
  VariantBuilder& set_integer_individual_field(const std::string& tag, const std::vector<int32_t>& values_for_all_samples);
  VariantBuilder& set_integer_individual_field(const std::string& tag, std::vector<int32_t>&& values_for_all_samples);
  VariantBuilder& set_integer_individual_field(const std::string& tag, const std::vector<std::vector<int32_t>>& values_for_all_samples);
  VariantBuilder& set_integer_individual_field(const std::string& tag, std::vector<std::vector<int32_t>>&& values_for_all_samples);
  VariantBuilder& set_integer_individual_field(const uint32_t field_index, const std::vector<int32_t>& values_for_all_samples);
  VariantBuilder& set_integer_individual_field(const uint32_t field_index, std::vector<int32_t>&& values_for_all_samples);
  VariantBuilder& set_integer_individual_field(const uint32_t field_index, const std::vector<std::vector<int32_t>>& values_for_all_samples);
  VariantBuilder& set_integer_individual_field(const uint32_t field_index, std::vector<std::vector<int32_t>>&& values_for_all_samples);
  VariantBuilder& set_float_individual_field(const std::string& tag, const std::vector<float>& values_for_all_samples);
  VariantBuilder& set_float_individual_field(const std::string& tag, std::vector<float>&& values_for_all_samples);
  VariantBuilder& set_float_individual_field(const std::string& tag, const std::vector<std::vector<float>>& values_for_all_samples);
  VariantBuilder& set_float_individual_field(const std::string& tag, std::vector<std::vector<float>>&& values_for_all_samples);
  VariantBuilder& set_float_individual_field(const uint32_t field_index, const std::vector<float>& values_for_all_samples);
  VariantBuilder& set_float_individual_field(const uint32_t field_index, std::vector<float>&& values_for_all_samples);
  VariantBuilder& set_float_individual_field(const uint32_t field_index, const std::vector<std::vector<float>>& values_for_all_samples);
  VariantBuilder& set_float_individual_field(const uint32_t field_index, std::vector<std::vector<float>>&& values_for_all_samples);
  VariantBuilder& set_string_individual_field(const std::string& tag, const std::vector<std::string>& values_for_all_samples);
  VariantBuilder& set_string_individual_field(const std::string& tag, std::vector<std::string>&& values_for_all_samples);
  VariantBuilder& set_string_individual_field(const uint32_t field_index, const std::vector<std::string>& values_for_all_samples);
  VariantBuilder& set_string_individual_field(const uint32_t field_index, std::vector<std::string>&& values_for_all_samples);

  /******************************************************************************
   *
   * Functions for setting individual/FORMAT fields by sample
   *
   * -These take values for one sample at a time, and so are less efficient
   *  than the bulk individual field setters above
   *
   * -Setting by field index and sample index is more efficient than setting
   *  by field name and sample name, provided that you look up the indices ONCE
   *  in the header and cache them instead of looking them up for every record.
   *
   * -The Genotype setter functions require you to encode your genotype data
   *  with one of the Genotype::encode_genotype() functions before passing it in.
   *
   *  For example, if you want to encode the genotype 0/1, create a vector with
   *  {0, 1}, pass it to Genotype::encode_genotype(), then pass it to the appropriate
   *  set_genotypes() function below.
   *
   ******************************************************************************/

  VariantBuilder& set_genotype(const std::string& sample, const std::vector<int32_t>& genotype);
  VariantBuilder& set_genotype(const uint32_t sample_index, const std::vector<int32_t>& genotype);
  VariantBuilder& set_integer_individual_field(const std::string& tag, const std::string& sample, const int32_t value);
  VariantBuilder& set_integer_individual_field(const std::string& tag, const std::string& sample, const std::vector<int32_t>& values);
  VariantBuilder& set_integer_individual_field(const uint32_t field_index, const uint32_t sample_index, const int32_t value);
  VariantBuilder& set_integer_individual_field(const uint32_t field_index, const uint32_t sample_index, const std::vector<int32_t>& values);
  VariantBuilder& set_float_individual_field(const std::string& tag, const std::string& sample, const float value);
  VariantBuilder& set_float_individual_field(const std::string& tag, const std::string& sample, const std::vector<float>& values);
  VariantBuilder& set_float_individual_field(const uint32_t field_index, const uint32_t sample_index, const float value);
  VariantBuilder& set_float_individual_field(const uint32_t field_index, const uint32_t sample_index, const std::vector<float>& values);
  VariantBuilder& set_string_individual_field(const std::string& tag, const std::string& sample, const std::string& value);
  VariantBuilder& set_string_individual_field(const uint32_t field_index, const uint32_t sample_index, const std::string& value);

  /******************************************************************************
   *
   * Functions for removing individual/FORMAT fields
   *
   ******************************************************************************/

  VariantBuilder& remove_individual_field(const std::string& tag);
  VariantBuilder& remove_individual_field(const uint32_t field_index);
  VariantBuilder& remove_individual_fields(const std::vector<std::string>& tags);
  VariantBuilder& remove_individual_fields(const std::vector<uint32_t>& field_indices);

  /******************************************************************************
   *
   * Builder operations: build() and clear()
   *
   * -In general you should create ONE VariantBuilder and re-use it across
   *  records, calling clear() in between records. This is much more efficient
   *  than creating a new builder for each record, since it greatly cuts down
   *  on memory allocations/de-allocations for the builder's internal data
   *  structures.
   *
   ******************************************************************************/

  Variant build() const;
  VariantBuilder& clear();

  /******************************************************************************
   *
   * Private functionality
   *
   ******************************************************************************/

 private:
  VariantHeader m_header;
  VariantBuilderCoreField<int32_t> m_contig;
  VariantBuilderCoreField<int32_t> m_start_pos;
  VariantBuilderCoreField<int32_t> m_stop_pos;
  VariantBuilderCoreField<float> m_qual;
  VariantBuilderSharedRegion m_shared_region;
  VariantBuilderIndividualRegion m_individual_region;
  bool m_enable_validation;

  void build_from_scratch(const std::shared_ptr<bcf1_t>& new_variant_body) const;
  void post_build_validation(const std::shared_ptr<bcf1_t>& new_variant_body) const;
};

}

#endif  /* gamgee__variant_builder__guard */
