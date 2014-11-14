#ifndef gamgee__variant_builder__guard
#define gamgee__variant_builder__guard

#include "variant.h"
#include "variant_builder_shared_region.h"
#include "variant_builder_individual_region.h"
#include "variant_builder_multi_sample_vector.h"

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
 * record(s) you intend to create (with all shared/individual fields,
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
 * See the comments further down for instructions and tips on using the various
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
 * ------------------------
 *
 * -Setting by field/sample index is more efficient than setting by field/sample name, provided
 *  that you look up the index for each field/sample ONCE in the header and cache it at traversal
 *  start instead of looking it up for every record.
 *
 * -For setting single-valued fields it's more efficient to use the functions that take a
 *  scalar value (int, float, etc.) instead of a vector.
 *
 * -For removing fields you generally have the option of either calling the appropriate remove_*
 *  API function, OR passing in a missing/empty value to the appropriate set_* function. Both
 *  options are equivalent and will result in the field being removed.
 *
 * -Disabling validation is possible, and will certainly improve performance, however if you
 *  attempt to perform an action that would have been prevented by validation checks (such as
 *  setting a non-existent field) you WILL get undefined behavior in your program. You should
 *  ONLY disable validation if you're extremely confident that your program will not take any
 *  invalid/incorrect actions, and that the data you pass to the builder will always be valid.
 *
 * Notes on setting individual fields:
 * -----------------------------------
 *
 * -Setting individual fields by move/r-value is more efficient than setting by l-value, as the
 *  functions that take an l-value reference assume that you want to keep your data and therefore
 *  make a copy of it. Note that you will usually have to invoke std::move() explicitly to avoid
 *  a copy -- eg.,
 *
 *      builder.set_integer_individual_field(field_index, std::move(my_vector));
 *
 * -Setting individual fields in bulk (ie., all samples at once) is more efficient than setting
 *  one sample at a time.
 *
 * -It is an error to request both bulk and per-sample changes to the same field
 *  (without calling clear() in between). This is because it would be too costly
 *  to reconcile the two kinds of changes.
 *
 * -The bulk-setting functions that take a VariantBuilderMultiSampleVector (which internally
 *  is a pre-padded one-dimensional vector of values) are much more efficient than the
 *  functions that take a two-dimensional vector, but the VariantBuilderMultiSampleVector
 *  approach requires a bit more work to use, and also requires that you know the maximum
 *  number of values per sample for the field in advance.
 *
 *  For example, with 4 samples and an integer individual field with a
 *  varied number of values per sample, you could pass in the following
 *  two-dimensional vector:
 *
 *  { {1, 2}, {3}, {}, {5, 6, 7} }
 *
 *  with each inner vector representing the values for one sample. However,
 *  this nested vector is a fairly inefficient data structure with poor
 *  data locality/cache performance.
 *
 *  If high performance is desired, you can use a VariantBuilderMultiSampleVector
 *  instead of a two-dimensional vector:
 *
 *  First determine the number of samples and the maximum number of values per
 *  sample for the field. Then get a pre-initialized VariantBuilderMultiSampleVector
 *  from the builder:
 *
 *  auto multi_sample_vector = builder.get_integer_multi_sample_vector(num_samples, max_values_per_sample);
 *
 *  This vector will have missing values for all samples, with appropriate padding to the maximum field width.
 *
 *  Then, fill in the values for each non-missing sample by invoking the set_sample_value() and/or
 *  set_sample_values() functions on your multi-sample vector (NOTE: set_sample_value() is MUCH more efficient
 *  than set_sample_values() since it doesn't require a vector construction/destruction for each call).
 *  You don't have to worry about samples with no values, since all samples start out with missing values.
 *
 *  Finally, pass your multi-sample vector into a setter function (favoring the functions
 *  that take field indices and use move semantics for high performance):
 *
 *  builder.set_integer_individual_field(field_index, std::move(multi_sample_vector));
 *
 *  The advantage of the VariantBuilderMultiSampleVector approach is much greater efficiency
 *  in terms of data locality and memory allocations.
 *
 * -The Genotype setter functions require you to encode your genotype data
 *  with one of the Genotype::encode_genotype() functions before passing it in.
 *
 *  For example, if you want to encode the genotype 0/1, create a vector with
 *  {0, 1}, pass it to Genotype::encode_genotype(), then pass it to the appropriate
 *  set_genotype() function.
 *
 *  Multi-sample collections of genotypes (either a VariantBuilderMultiSampleVector or
 *  a vector<vector>) must also be encoded before passing them to the builder. Use
 *  Genotype::encode_genotypes() to encode them (note the plural).
 *
 */
class VariantBuilder {
 public:

  /**
   * @brief Construct a new VariantBuilder given a VariantHeader
   *
   * @param header VariantHeader to use for constructing and validating new Variant records
   *
   * @note header will become the header for all Variant records created, and so must be appropriate for the final
   *       Variant records you intend to create (with all shared/individual fields, samples, and contigs pre-declared)
   *
   * @note It's much more efficient to create a single VariantBuilder and re-use it across records (calling clear()
   *       as needed) instead of creating a new VariantBuilder for each record. VariantBuilder construction involves
   *       many memory allocations and is not cheap!
   */
  explicit VariantBuilder(const VariantHeader& header);

  // VariantBuilders are moveable but not copyable, and use default destruction
  VariantBuilder(VariantBuilder&& other) = default;
  VariantBuilder& operator=(VariantBuilder&& other) = default;
  VariantBuilder(const VariantBuilder& other) = delete;
  VariantBuilder& operator=(const VariantBuilder& other) = delete;
  ~VariantBuilder() = default;

  /**
   * @brief Return the Variant header for this builder
   *
   * @return Variant header this builder was constructed with (and which becomes the header for any Variant records built)
   */
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

  /**
   * @brief Set the chromosome by index
   *
   * @param chromosome chromosome index (from a header lookup)
   */
  VariantBuilder& set_chromosome(const uint32_t chromosome);

  /**
   * @brief Set the chromosome by name
   *
   * @param chromosome chromosome name
   */
  VariantBuilder& set_chromosome(const std::string& chromosome);

  /**
   * @brief Set the alignment start position
   *
   * @param alignment_start 1-based alignment start position (as you would see in a VCF file)
   *
   * @note the internal encoding is 0-based to mimic that of the BCF files
   */
  VariantBuilder& set_alignment_start(const uint32_t alignment_start);

  /**
   * @brief Set the alignment stop position
   *
   * @param alignment_stop 1-based alignment stop position (as you would see in a VCF INFO END tag)
   *
   * @note the internal encoding is 0-based to mimic that of the BCF files
   */
  VariantBuilder& set_alignment_stop(const uint32_t alignment_stop);

  /**
   * @brief Set the Phred-scaled site quality (probability that the site is not reference)
   *
   * @param qual Phred-scaled site quality (probability that the site is not reference)
   */
  VariantBuilder& set_qual(const float qual);

  /**
   * @brief Set the variant ID field
   *
   * @param id variant ID (typically DBSNP ID)
   */
  VariantBuilder& set_id(const std::string& id);

  /**
   * @brief Set the reference allele
   *
   * @param ref_allele reference allele as a string
   */
  VariantBuilder& set_ref_allele(const std::string& ref_allele);

  /**
   * @brief Set the alt allele
   *
   * @param alt_allele alt allele as a string
   *
   * @note Alt allele passed in replaces any previous value(s) for the alt field
   * @note It's more efficient to use this function when there's only a single alt allele
   *       instead of the vector-based setter
   */
  VariantBuilder& set_alt_allele(const std::string& alt_allele);

  /**
   * @brief Set the alt alleles
   *
   * @param alt_alleles one string per alt allele
   *
   * @note Alt alleles passed in replace any previous value(s) for the alt field
   * @note It's more efficient to use the setter that takes a single string when there's
   *       only one alt allele
   */
  VariantBuilder& set_alt_alleles(const std::vector<std::string>& alt_alleles);

  /**
   * @brief Set the filters using filter names
   *
   * @param filters vector of filter names
   *
   * @note It's more efficient to set using filter indices
   */
  VariantBuilder& set_filters(const std::vector<std::string>& filters);

  /**
   * @brief Set the filters using filter indices
   *
   * @param filters vector of filter indices (from header lookups)
   */
  VariantBuilder& set_filters(const std::vector<int32_t>& filters);


  /******************************************************************************
   *
   * Functions for removing core site-level fields
   * (Note that we don't have removal functions for the required fields)
   *
   ******************************************************************************/

  /**
   * @brief Clear the alignment stop value (if set)
   */
  VariantBuilder& remove_alignment_stop();

  /**
   * @brief Clear the qual field (if set)
   */
  VariantBuilder& remove_qual();

  /**
   * @brief Clear the ID field (if set)
   */
  VariantBuilder& remove_id();

  /**
   * @brief Clear alt alleles (if set)
   */
  VariantBuilder& remove_alt_alleles();

  /**
   * @brief Clear filters (if set)
   */
  VariantBuilder& remove_filters();


  /******************************************************************************
   *
   * Functions for setting shared/INFO fields
   *
   ******************************************************************************/

  /**
   * @brief Set a single-valued integer shared field by field name
   *
   * @param tag name of the shared field to set
   * @param value new value for the field
   *
   * @note It's more efficient to use this setter instead of a vector-based setter when
   *       a field has just a single value
   *
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_integer_shared_field(const std::string& tag, const int32_t value);

  /**
   * @brief Set a multi-valued integer shared field by field name
   *
   * @param tag name of the shared field to set
   * @param values new values for the field
   *
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_integer_shared_field(const std::string& tag, const std::vector<int32_t>& values);

  /**
   * @brief Set a single-valued integer shared field by field index
   *
   * @param index index of the shared field to set (from a header lookup)
   * @param value new value for the field
   *
   * @note It's more efficient to use this setter instead of a vector-based setter when
   *       a field has just a single value
   */
  VariantBuilder& set_integer_shared_field(const uint32_t index, const int32_t value);

  /**
   * @brief Set a multi-valued integer shared field by field index
   *
   * @param index index of the shared field to set (from a header lookup)
   * @param values new values for the field
   */
  VariantBuilder& set_integer_shared_field(const uint32_t index, const std::vector<int32_t>& values);

  /**
   * @brief Set a single-valued float shared field by field name
   *
   * @param tag name of the shared field to set
   * @param value new value for the field
   *
   * @note It's more efficient to use this setter instead of a vector-based setter when
   *       a field has just a single value
   *
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_float_shared_field(const std::string& tag, const float value);

  /**
   * @brief Set a multi-valued float shared field by field name
   *
   * @param tag name of the shared field to set
   * @param values new values for the field
   *
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_float_shared_field(const std::string& tag, const std::vector<float>& values);

  /**
   * @brief Set a single-valued float shared field by field index
   *
   * @param index index of the shared field to set (from a header lookup)
   * @param value new value for the field
   *
   * @note It's more efficient to use this setter instead of a vector-based setter when
   *       a field has just a single value
   */
  VariantBuilder& set_float_shared_field(const uint32_t index, const float value);

  /**
   * @brief Set a multi-valued float shared field by field index
   *
   * @param index index of the shared field to set (from a header lookup)
   * @param values new values for the field
   */
  VariantBuilder& set_float_shared_field(const uint32_t index, const std::vector<float>& values);

  /**
   * @brief Set a string shared field by field name
   *
   * @param tag name of the shared field to set
   * @param value new value for the field
   *
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_string_shared_field(const std::string& tag, const std::string& value);

  /**
   * @brief Set a string shared field by field index
   *
   * @param index index of the shared field to set (from a header lookup)
   * @param value new value for the field
   */
  VariantBuilder& set_string_shared_field(const uint32_t index, const std::string& value);

  /**
   * @brief Set a boolean (flag) shared field by field name
   *
   * @param tag name of the shared field to set
   *
   * @note Field is set to true/present. To set a boolean field to false call one of the
   *       remove_shared_field() functions.
   *
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_boolean_shared_field(const std::string& tag);

  /**
   * @brief Set a boolean (flag) shared field by field index
   *
   * @param index index of the shared field to set (from a header lookup)
   *
   * @note Field is set to true/present. To set a boolean field to false call one of the
   *       remove_shared_field() functions.
   */
  VariantBuilder& set_boolean_shared_field(const uint32_t index);


  /******************************************************************************
   *
   * Functions for removing shared/INFO fields
   *
   ******************************************************************************/

  /**
   * @brief Remove a shared field by field name
   *
   * @param tag name of the shared field to remove
   *
   * @note Less efficient than removing using the field index
   */
  VariantBuilder& remove_shared_field(const std::string& tag);

  /**
   * @brief Remove a shared field by field index
   *
   * @param field_index index of the shared field to remove (from a header lookup)
   */
  VariantBuilder& remove_shared_field(const uint32_t field_index);

  /**
   * @brief Remove multiple shared fields by field name
   *
   * @param tags names of the shared fields to remove
   *
   * @note Less efficient than removing using the field indices
   */
  VariantBuilder& remove_shared_fields(const std::vector<std::string>& tags);

  /**
   * @brief Remove multiple shared fields by field index
   *
   * @param field_indices indices of the shared fields to remove (from header lookups)
   */
  VariantBuilder& remove_shared_fields(const std::vector<uint32_t>& field_indices);


  /******************************************************************************
   *
   * Functions for setting individual/FORMAT fields in bulk
   * (ie., setting all values for all samples at once)
   *
   ******************************************************************************/

  /**
   * @brief Set the genotypes (GT) field for all samples at once using an efficient flattened (one-dimensional) vector, COPYING the provided values
   *
   * @param genotypes_for_all_samples encoded genotypes for all samples as a VariantBuilderMultiSampleVector (see note below)
   *
   * @note must pass vector to appropriate Genotype::encode_genotypes() function before passing it to this function
   *
   * @note To create a multi-sample flattened vector for use with this function, first determine the number of samples
   *       and the maximum ploidy across all samples, then get a pre-initialized vector from the builder:
   *
   *       auto multi_sample_vector = builder.get_genotype_multi_sample_vector(num_samples, max_ploidy);
   *
   *       This vector will have missing values for all samples, with appropriate padding to the maximum ploidy.
   *
   *       Then, fill in the values for each non-missing sample by invoking the set_sample_value() and/or
   *       set_sample_values() functions on your multi-sample vector (set_sample_value() is more efficient
   *       than set_sample_values() since it doesn't require a vector construction/destruction for each call).
   *       You don't have to worry about samples with no values, since all samples start out with missing values,
   *       however you should represent each no-call allele with -1.
   *
   *       Finally, pass your multi-sample vector into this function:
   *
   *       builder.set_genotypes(multi_sample_vector);
   *
   *       If this process is too inconvenient, or you can't know the maximum ploidy in advance,
   *       use a less-efficient function that takes a nested vector.
   *
   * @note Less efficient than moving the values into the builder
   */
  VariantBuilder& set_genotypes(const VariantBuilderMultiSampleVector<int32_t>& genotypes_for_all_samples);

  /**
   * @brief Set the genotypes (GT) field for all samples at once using an efficient flattened (one-dimensional) vector, MOVING the provided values
   *
   * @param genotypes_for_all_samples encoded genotypes for all samples as a VariantBuilderMultiSampleVector (see note below)
   *
   * @note must pass vector to appropriate Genotype::encode_genotypes() function before passing it to this function
   *
   * @note To create a multi-sample flattened vector for use with this function, first determine the number of samples
   *       and the maximum ploidy across all samples, then get a pre-initialized vector from the builder:
   *
   *       auto multi_sample_vector = builder.get_genotype_multi_sample_vector(num_samples, max_ploidy);
   *
   *       This vector will have missing values for all samples, with appropriate padding to the maximum ploidy.
   *
   *       Then, fill in the values for each non-missing sample by invoking the set_sample_value() and/or
   *       set_sample_values() functions on your multi-sample vector (set_sample_value() is more efficient
   *       than set_sample_values() since it doesn't require a vector construction/destruction for each call).
   *       You don't have to worry about samples with no values, since all samples start out with missing values,
   *       however you should represent each no-call allele with -1.
   *
   *       Finally, MOVE your multi-sample vector into this function:
   *
   *       builder.set_genotypes(std::move(multi_sample_vector));
   *
   *       If this process is too inconvenient, or you can't know the maximum ploidy in advance,
   *       use a less-efficient function that takes a nested vector.
   */
  VariantBuilder& set_genotypes(VariantBuilderMultiSampleVector<int32_t>&& genotypes_for_all_samples);

  /**
   * @brief Set the genotypes (GT) field for all samples at once by nested vector, copying the provided values
   *
   * @param genotypes_for_all_samples encoded genotypes for all samples in order of sample index, with one inner vector per sample (no padding necessary)
   *
   * @note must pass vector to appropriate Genotype::encode_genotypes() function before passing it to this function
   *
   * @note With a nested vector, each inner vector represents the genotypes for the sample with the corresponding index.
   *       There is no need to manually pad to the maximum ploidy, but you do need to add a missing value (-1) for
   *       each missing/no-call allele.
   *
   *       For example, if you had Sample1=0/1 Sample2=./. Sample3=. Sample4=0/1/2
   *       you would need to create the following nested vector, then pass it to Genotype::encode_genotypes():
   *
   *       { {0, 1}, {-1, -1}, {}, {0, 1, 2} }
   *
   * @note Less efficient than using a flattened vector
   * @note Less efficient than moving the values into the builder
   */
  VariantBuilder& set_genotypes(const std::vector<std::vector<int32_t>>& genotypes_for_all_samples);

  /**
   * @brief Set the genotypes (GT) field for all samples at once by nested vector, moving the provided values
   *
   * @param genotypes_for_all_samples encoded genotypes for all samples in order of sample index, with one inner vector per sample (no padding necessary)
   *
   * @note must pass vector to appropriate Genotype::encode_genotypes() function before passing it to this function
   *
   * @note With a nested vector, each inner vector represents the genotypes for the sample with the corresponding index.
   *       There is no need to manually pad to the maximum ploidy, but you do need to add a missing value (-1) for
   *       each missing/no-call allele.
   *
   *       For example, if you had Sample1=0/1 Sample2=./. Sample3=. Sample4=0/1/2
   *       you would need to create the following nested vector, then pass it to Genotype::encode_genotypes():
   *
   *       { {0, 1}, {-1, -1}, {}, {0, 1, 2} }
   *
   * @note Less efficient than using a flattened vector
   */
  VariantBuilder& set_genotypes(std::vector<std::vector<int32_t>>&& genotypes_for_all_samples);

  /**
   * @brief Set an integer individual field for all samples at once by name using an efficient flattened (one-dimensional) vector, COPYING the provided values
   *
   * @param tag name of the individual field to set
   * @param values_for_all_samples field values for all samples as a VariantBuilderMultiSampleVector (see note below)
   *
   * @note To create a multi-sample flattened vector for use with this function, first determine the number of samples
   *       and the maximum number of values per sample for this field, then get a pre-initialized vector from the builder:
   *
   *       auto multi_sample_vector = builder.get_integer_multi_sample_vector(num_samples, max_values_per_sample);
   *
   *       This vector will have missing values for all samples, with appropriate padding to the maximum field width.
   *
   *       Then, fill in the values for each non-missing sample by invoking the set_sample_value() and/or
   *       set_sample_values() functions on your multi-sample vector (set_sample_value() is more efficient
   *       than set_sample_values() since it doesn't require a vector construction/destruction for each call).
   *       You don't have to worry about samples with no values, since all samples start out with missing values.
   *
   *       Finally, pass your multi-sample vector into this function:
   *
   *       builder.set_integer_individual_field(field_name, multi_sample_vector);
   *
   *       If this process is too inconvenient, or you can't know the maximum number of values per sample in advance,
   *       use a less-efficient function that takes a nested vector.
   *
   * @note Less efficient than moving the values into the builder
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_integer_individual_field(const std::string& tag, const VariantBuilderMultiSampleVector<int32_t>& values_for_all_samples);

  /**
   * @brief Set an integer individual field for all samples at once by name using an efficient flattened (one-dimensional) vector, MOVING the provided values
   *
   * @param tag name of the individual field to set
   * @param values_for_all_samples field values for all samples as a VariantBuilderMultiSampleVector (see note below)
   *
   * @note To create a multi-sample flattened vector for use with this function, first determine the number of samples
   *       and the maximum number of values per sample for this field, then get a pre-initialized vector from the builder:
   *
   *       auto multi_sample_vector = builder.get_integer_multi_sample_vector(num_samples, max_values_per_sample);
   *
   *       This vector will have missing values for all samples, with appropriate padding to the maximum field width.
   *
   *       Then, fill in the values for each non-missing sample by invoking the set_sample_value() and/or
   *       set_sample_values() functions on your multi-sample vector (set_sample_value() is more efficient
   *       than set_sample_values() since it doesn't require a vector construction/destruction for each call).
   *       You don't have to worry about samples with no values, since all samples start out with missing values.
   *
   *       Finally, MOVE your multi-sample vector into this function:
   *
   *       builder.set_integer_individual_field(field_name, std::move(multi_sample_vector));
   *
   *       If this process is too inconvenient, or you can't know the maximum number of values per sample in advance,
   *       use a less-efficient function that takes a nested vector.
   *
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_integer_individual_field(const std::string& tag, VariantBuilderMultiSampleVector<int32_t>&& values_for_all_samples);

  /**
   * @brief Set an integer individual field for all samples at once by name using a nested vector, copying the provided values
   *
   * @param tag name of the individual field to set
   * @param values_for_all_samples field values for all samples in order of sample index, with one inner vector per sample (no special padding necessary)
   *
   * @note With a nested vector, each inner vector represents the values for the sample with the corresponding index.
   *       There is no need to manually pad with missing/vector end values. For example:
   *
   *       { {1, 2}, {3}, {}, {5, 6, 7} }
   *
   * @note Less efficient than using a flattened vector
   * @note Less efficient than moving the values into the builder
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_integer_individual_field(const std::string& tag, const std::vector<std::vector<int32_t>>& values_for_all_samples);

  /**
   * @brief Set an integer individual field for all samples at once by name using a nested vector, moving the provided values
   *
   * @param tag name of the individual field to set
   * @param values_for_all_samples field values for all samples in order of sample index, with one inner vector per sample (no special padding necessary)
   *
   * @note With a nested vector, each inner vector represents the values for the sample with the corresponding index.
   *       There is no need to manually pad with missing/vector end values. For example:
   *
   *       { {1, 2}, {3}, {}, {5, 6, 7} }
   *
   * @note Less efficient than using a flattened vector
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_integer_individual_field(const std::string& tag, std::vector<std::vector<int32_t>>&& values_for_all_samples);

  /**
   * @brief Set an integer individual field for all samples at once by index using an efficient flattened (one-dimensional) vector, COPYING the provided values
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param values_for_all_samples field values for all samples as a VariantBuilderMultiSampleVector (see note below)
   *
   * @note To create a multi-sample flattened vector for use with this function, first determine the number of samples
   *       and the maximum number of values per sample for this field, then get a pre-initialized vector from the builder:
   *
   *       auto multi_sample_vector = builder.get_integer_multi_sample_vector(num_samples, max_values_per_sample);
   *
   *       This vector will have missing values for all samples, with appropriate padding to the maximum field width.
   *
   *       Then, fill in the values for each non-missing sample by invoking the set_sample_value() and/or
   *       set_sample_values() functions on your multi-sample vector (set_sample_value() is more efficient
   *       than set_sample_values() since it doesn't require a vector construction/destruction for each call).
   *       You don't have to worry about samples with no values, since all samples start out with missing values.
   *
   *       Finally, pass your multi-sample vector into this function:
   *
   *       builder.set_integer_individual_field(field_index, multi_sample_vector);
   *
   *       If this process is too inconvenient, or you can't know the maximum number of values per sample in advance,
   *       use a less-efficient function that takes a nested vector.
   *
   * @note Less efficient than moving the values into the builder
   */
  VariantBuilder& set_integer_individual_field(const uint32_t field_index, const VariantBuilderMultiSampleVector<int32_t>& values_for_all_samples);

  /**
   * @brief Set an integer individual field for all samples at once by index using an efficient flattened (one-dimensional) vector, MOVING the provided values
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param values_for_all_samples field values for all samples as a VariantBuilderMultiSampleVector (see note below)
   *
   * @note To create a multi-sample flattened vector for use with this function, first determine the number of samples
   *       and the maximum number of values per sample for this field, then get a pre-initialized vector from the builder:
   *
   *       auto multi_sample_vector = builder.get_integer_multi_sample_vector(num_samples, max_values_per_sample);
   *
   *       This vector will have missing values for all samples, with appropriate padding to the maximum field width.
   *
   *       Then, fill in the values for each non-missing sample by invoking the set_sample_value() and/or
   *       set_sample_values() functions on your multi-sample vector (set_sample_value() is more efficient
   *       than set_sample_values() since it doesn't require a vector construction/destruction for each call).
   *       You don't have to worry about samples with no values, since all samples start out with missing values.
   *
   *       Finally, MOVE your multi-sample vector into this function:
   *
   *       builder.set_integer_individual_field(field_index, std::move(multi_sample_vector));
   *
   *       If this process is too inconvenient, or you can't know the maximum number of values per sample in advance,
   *       use a less-efficient function that takes a nested vector.
   */
  VariantBuilder& set_integer_individual_field(const uint32_t field_index, VariantBuilderMultiSampleVector<int32_t>&& values_for_all_samples);

  /**
   * @brief Set an integer individual field for all samples at once by index using a nested vector, copying the provided values
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param values_for_all_samples field values for all samples in order of sample index, with one inner vector per sample (no special padding necessary)
   *
   * @note With a nested vector, each inner vector represents the values for the sample with the corresponding index.
   *       There is no need to manually pad with missing/vector end values. For example:
   *
   *       { {1, 2}, {3}, {}, {5, 6, 7} }
   *
   * @note Less efficient than using a flattened vector
   * @note Less efficient than moving the values into the builder
   */
  VariantBuilder& set_integer_individual_field(const uint32_t field_index, const std::vector<std::vector<int32_t>>& values_for_all_samples);

  /**
   * @brief Set an integer individual field for all samples at once by index using a nested vector, moving the provided values
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param values_for_all_samples field values for all samples in order of sample index, with one inner vector per sample (no special padding necessary)
   *
   * @note With a nested vector, each inner vector represents the values for the sample with the corresponding index.
   *       There is no need to manually pad with missing/vector end values. For example:
   *
   *       { {1, 2}, {3}, {}, {5, 6, 7} }
   *
   * @note Less efficient than using a flattened vector
   */
  VariantBuilder& set_integer_individual_field(const uint32_t field_index, std::vector<std::vector<int32_t>>&& values_for_all_samples);

  /**
   * @brief Set a float individual field for all samples at once by name using an efficient flattened (one-dimensional) vector, COPYING the provided values
   *
   * @param tag name of the individual field to set
   * @param values_for_all_samples field values for all samples as a VariantBuilderMultiSampleVector (see note below)
   *
   * @note To create a multi-sample flattened vector for use with this function, first determine the number of samples
   *       and the maximum number of values per sample for this field, then get a pre-initialized vector from the builder:
   *
   *       auto multi_sample_vector = builder.get_float_multi_sample_vector(num_samples, max_values_per_sample);
   *
   *       This vector will have missing values for all samples, with appropriate padding to the maximum field width.
   *
   *       Then, fill in the values for each non-missing sample by invoking the set_sample_value() and/or
   *       set_sample_values() functions on your multi-sample vector (set_sample_value() is more efficient
   *       than set_sample_values() since it doesn't require a vector construction/destruction for each call).
   *       You don't have to worry about samples with no values, since all samples start out with missing values.
   *
   *       Finally, pass your multi-sample vector into this function:
   *
   *       builder.set_float_individual_field(field_name, multi_sample_vector);
   *
   *       If this process is too inconvenient, or you can't know the maximum number of values per sample in advance,
   *       use a less-efficient function that takes a nested vector.
   *
   * @note Less efficient than moving the values into the builder
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_float_individual_field(const std::string& tag, const VariantBuilderMultiSampleVector<float>& values_for_all_samples);

  /**
   * @brief Set a float individual field for all samples at once by name using an efficient flattened (one-dimensional) vector, MOVING the provided values
   *
   * @param tag name of the individual field to set
   * @param values_for_all_samples field values for all samples as a VariantBuilderMultiSampleVector (see note below)
   *
   * @note To create a multi-sample flattened vector for use with this function, first determine the number of samples
   *       and the maximum number of values per sample for this field, then get a pre-initialized vector from the builder:
   *
   *       auto multi_sample_vector = builder.get_float_multi_sample_vector(num_samples, max_values_per_sample);
   *
   *       This vector will have missing values for all samples, with appropriate padding to the maximum field width.
   *
   *       Then, fill in the values for each non-missing sample by invoking the set_sample_value() and/or
   *       set_sample_values() functions on your multi-sample vector (set_sample_value() is more efficient
   *       than set_sample_values() since it doesn't require a vector construction/destruction for each call).
   *       You don't have to worry about samples with no values, since all samples start out with missing values.
   *
   *       Finally, MOVE your multi-sample vector into this function:
   *
   *       builder.set_float_individual_field(field_name, std::move(multi_sample_vector));
   *
   *       If this process is too inconvenient, or you can't know the maximum number of values per sample in advance,
   *       use a less-efficient function that takes a nested vector.
   *
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_float_individual_field(const std::string& tag, VariantBuilderMultiSampleVector<float>&& values_for_all_samples);

  /**
   * @brief Set a float individual field for all samples at once by name using a nested vector, copying the provided values
   *
   * @param tag name of the individual field to set
   * @param values_for_all_samples field values for all samples in order of sample index, with one inner vector per sample (no special padding necessary)
   *
   * @note With a nested vector, each inner vector represents the values for the sample with the corresponding index.
   *       There is no need to manually pad with missing/vector end values. For example:
   *
   *       { {1.5, 2.5}, {3.5}, {}, {5.5, 6.5, 7.5} }
   *
   * @note Less efficient than using a flattened vector
   * @note Less efficient than moving the values into the builder
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_float_individual_field(const std::string& tag, const std::vector<std::vector<float>>& values_for_all_samples);

  /**
   * @brief Set a float individual field for all samples at once by name using a nested vector, moving the provided values
   *
   * @param tag name of the individual field to set
   * @param values_for_all_samples field values for all samples in order of sample index, with one inner vector per sample (no special padding necessary)
   *
   * @note With a nested vector, each inner vector represents the values for the sample with the corresponding index.
   *       There is no need to manually pad with missing/vector end values. For example:
   *
   *       { {1.5, 2.5}, {3.5}, {}, {5.5, 6.5, 7.5} }
   *
   * @note Less efficient than using a flattened vector
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_float_individual_field(const std::string& tag, std::vector<std::vector<float>>&& values_for_all_samples);

  /**
   * @brief Set a float individual field for all samples at once by index using an efficient flattened (one-dimensional) vector, COPYING the provided values
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param values_for_all_samples field values for all samples as a VariantBuilderMultiSampleVector (see note below)
   *
   * @note To create a multi-sample flattened vector for use with this function, first determine the number of samples
   *       and the maximum number of values per sample for this field, then get a pre-initialized vector from the builder:
   *
   *       auto multi_sample_vector = builder.get_float_multi_sample_vector(num_samples, max_values_per_sample);
   *
   *       This vector will have missing values for all samples, with appropriate padding to the maximum field width.
   *
   *       Then, fill in the values for each non-missing sample by invoking the set_sample_value() and/or
   *       set_sample_values() functions on your multi-sample vector (set_sample_value() is more efficient
   *       than set_sample_values() since it doesn't require a vector construction/destruction for each call).
   *       You don't have to worry about samples with no values, since all samples start out with missing values.
   *
   *       Finally, pass your multi-sample vector into this function:
   *
   *       builder.set_float_individual_field(field_index, multi_sample_vector);
   *
   *       If this process is too inconvenient, or you can't know the maximum number of values per sample in advance,
   *       use a less-efficient function that takes a nested vector.
   *
   * @note Less efficient than moving the values into the builder
   */
  VariantBuilder& set_float_individual_field(const uint32_t field_index, const VariantBuilderMultiSampleVector<float>& values_for_all_samples);

  /**
   * @brief Set a float individual field for all samples at once by index using an efficient flattened (one-dimensional) vector, MOVING the provided values
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param values_for_all_samples field values for all samples as a VariantBuilderMultiSampleVector (see note below)
   *
   * @note To create a multi-sample flattened vector for use with this function, first determine the number of samples
   *       and the maximum number of values per sample for this field, then get a pre-initialized vector from the builder:
   *
   *       auto multi_sample_vector = builder.get_float_multi_sample_vector(num_samples, max_values_per_sample);
   *
   *       This vector will have missing values for all samples, with appropriate padding to the maximum field width.
   *
   *       Then, fill in the values for each non-missing sample by invoking the set_sample_value() and/or
   *       set_sample_values() functions on your multi-sample vector (set_sample_value() is more efficient
   *       than set_sample_values() since it doesn't require a vector construction/destruction for each call).
   *       You don't have to worry about samples with no values, since all samples start out with missing values.
   *
   *       Finally, MOVE your multi-sample vector into this function:
   *
   *       builder.set_float_individual_field(field_index, std::move(multi_sample_vector));
   *
   *       If this process is too inconvenient, or you can't know the maximum number of values per sample in advance,
   *       use a less-efficient function that takes a nested vector.
   */
  VariantBuilder& set_float_individual_field(const uint32_t field_index, VariantBuilderMultiSampleVector<float>&& values_for_all_samples);

  /**
   * @brief Set a float individual field for all samples at once by index using a nested vector, copying the provided values
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param values_for_all_samples field values for all samples in order of sample index, with one inner vector per sample (no special padding necessary)
   *
   * @note With a nested vector, each inner vector represents the values for the sample with the corresponding index.
   *       There is no need to manually pad with missing/vector end values. For example:
   *
   *       { {1.5, 2.5}, {3.5}, {}, {5.5, 6.5, 7.5} }
   *
   * @note Less efficient than using a flattened vector
   * @note Less efficient than moving the values into the builder
   */
  VariantBuilder& set_float_individual_field(const uint32_t field_index, const std::vector<std::vector<float>>& values_for_all_samples);

  /**
   * @brief Set a float individual field for all samples at once by index using a nested vector, moving the provided values
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param values_for_all_samples field values for all samples in order of sample index, with one inner vector per sample (no special padding necessary)
   *
   * @note With a nested vector, each inner vector represents the values for the sample with the corresponding index.
   *       There is no need to manually pad with missing/vector end values. For example:
   *
   *       { {1.5, 2.5}, {3.5}, {}, {5.5, 6.5, 7.5} }
   *
   * @note Less efficient than using a flattened vector
   */
  VariantBuilder& set_float_individual_field(const uint32_t field_index, std::vector<std::vector<float>>&& values_for_all_samples);

  /**
   * @brief Set a string individual field for all samples at once by name, copying the provided values
   *
   * @param tag name of the individual field to set
   * @param values_for_all_samples field values for all samples in order of sample index
   *
   * @note Less efficient than moving the values into the builder
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_string_individual_field(const std::string& tag, const std::vector<std::string>& values_for_all_samples);

  /**
   * @brief Set a string individual field for all samples at once by name, moving the provided values
   *
   * @param tag name of the individual field to set
   * @param values_for_all_samples field values for all samples in order of sample index
   *
   * @note Less efficient than setting using the field index
   */
  VariantBuilder& set_string_individual_field(const std::string& tag, std::vector<std::string>&& values_for_all_samples);

  /**
   * @brief Set a string individual field for all samples at once by index, copying the provided values
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param values_for_all_samples field values for all samples in order of sample index
   *
   * @note Less efficient than moving the values into the builder
   */
  VariantBuilder& set_string_individual_field(const uint32_t field_index, const std::vector<std::string>& values_for_all_samples);

  /**
   * @brief Set a string individual field for all samples at once by index, moving the provided values
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param values_for_all_samples field values for all samples in order of sample index
   */
  VariantBuilder& set_string_individual_field(const uint32_t field_index, std::vector<std::string>&& values_for_all_samples);


  /******************************************************************************
   *
   * Functions for setting individual/FORMAT fields by sample
   *
   ******************************************************************************/

  /**
   * @brief Set the genotypes (GT) field for a single sample by sample name
   *
   * @param sample name of the sample whose genotype to set
   * @param genotype encoded genotype for the specified sample (see notes below)
   *
   * Examples: For genotype 0/1, create vector {0, 1}, then pass it to Genotype::encode_genotype()
   *           For genotype 1/., create vector {1, -1}, then pass it to Genotype::encode_genotype()
   *           For genotype ./., create vector {-1, -1}, then pass it to Genotype::encode_genotype()
   *
   * @note must pass vector to appropriate Genotype::encode_genotype() function before passing it to this function
   *
   * @note Less efficient than setting using the sample index
   */
  VariantBuilder& set_genotype(const std::string& sample, const std::vector<int32_t>& genotype);

  /**
   * @brief Set the genotypes (GT) field for a single sample by sample index
   *
   * @param sample_index index of the sample whose genotype to set (from a header lookup)
   * @param genotype encoded genotype for the specified sample (see notes below)
   *
   * Examples: For genotype 0/1, create vector {0, 1}, then pass it to Genotype::encode_genotype()
   *           For genotype 1/., create vector {1, -1}, then pass it to Genotype::encode_genotype()
   *           For genotype ./., create vector {-1, -1}, then pass it to Genotype::encode_genotype()
   *
   * @note must pass vector to appropriate Genotype::encode_genotype() function before passing it to this function
   */
  VariantBuilder& set_genotype(const uint32_t sample_index, const std::vector<int32_t>& genotype);

  /**
   * @brief Set a single-valued integer individual field for a single sample by field and sample name
   *
   * @param tag name of the individual field to set
   * @param sample name of the sample whose value to set
   * @param value field value for the specified sample
   *
   * @note It's more efficient to use this setter instead of a vector-based setter when
   *       a sample has just a single value for a field
   * @note Less efficient than setting using the field/sample indices
   */
  VariantBuilder& set_integer_individual_field(const std::string& tag, const std::string& sample, const int32_t value);

  /**
   * @brief Set a multi-valued integer individual field for a single sample by field and sample name
   *
   * @param tag name of the individual field to set
   * @param sample name of the sample whose value to set
   * @param values field values for the specified sample
   *
   * @note Less efficient than setting using the field/sample indices
   */
  VariantBuilder& set_integer_individual_field(const std::string& tag, const std::string& sample, const std::vector<int32_t>& values);

  /**
   * @brief Set a single-valued integer individual field for a single sample by field and sample index
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param sample_index index of the sample whose value to set (from a header lookup)
   * @param value field value for the specified sample
   *
   * @note It's more efficient to use this setter instead of a vector-based setter when
   *       a sample has just a single value for a field
   */
  VariantBuilder& set_integer_individual_field(const uint32_t field_index, const uint32_t sample_index, const int32_t value);

  /**
   * @brief Set a multi-valued integer individual field for a single sample by field and sample index
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param sample_index index of the sample whose value to set (from a header lookup)
   * @param values field values for the specified sample
   */
  VariantBuilder& set_integer_individual_field(const uint32_t field_index, const uint32_t sample_index, const std::vector<int32_t>& values);

  /**
   * @brief Set a single-valued float individual field for a single sample by field and sample name
   *
   * @param tag name of the individual field to set
   * @param sample name of the sample whose value to set
   * @param value field value for the specified sample
   *
   * @note It's more efficient to use this setter instead of a vector-based setter when
   *       a sample has just a single value for a field
   * @note Less efficient than setting using the field/sample indices
   */
  VariantBuilder& set_float_individual_field(const std::string& tag, const std::string& sample, const float value);

  /**
   * @brief Set a multi-valued float individual field for a single sample by field and sample name
   *
   * @param tag name of the individual field to set
   * @param sample name of the sample whose value to set
   * @param values field values for the specified sample
   *
   * @note Less efficient than setting using the field/sample indices
   */
  VariantBuilder& set_float_individual_field(const std::string& tag, const std::string& sample, const std::vector<float>& values);

  /**
   * @brief Set a single-valued float individual field for a single sample by field and sample index
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param sample_index index of the sample whose value to set (from a header lookup)
   * @param value field value for the specified sample
   *
   * @note It's more efficient to use this setter instead of a vector-based setter when
   *       a sample has just a single value for a field
   */
  VariantBuilder& set_float_individual_field(const uint32_t field_index, const uint32_t sample_index, const float value);

  /**
   * @brief Set a multi-valued float individual field for a single sample by field and sample index
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param sample_index index of the sample whose value to set (from a header lookup)
   * @param values field values for the specified sample
   */
  VariantBuilder& set_float_individual_field(const uint32_t field_index, const uint32_t sample_index, const std::vector<float>& values);

  /**
   * @brief Set a string individual field for a single sample by field and sample name
   *
   * @param tag name of the individual field to set
   * @param sample name of the sample whose value to set
   * @param value field value for the specified sample
   *
   * @note Less efficient than setting using the field/sample indices
   */
  VariantBuilder& set_string_individual_field(const std::string& tag, const std::string& sample, const std::string& value);

  /**
   * @brief Set a string individual field for a single sample by field and sample index
   *
   * @param field_index index of the individual field to set (from a header lookup)
   * @param sample_index index of the sample whose value to set (from a header lookup)
   * @param value field value for the specified sample
   */
  VariantBuilder& set_string_individual_field(const uint32_t field_index, const uint32_t sample_index, const std::string& value);


  /******************************************************************************
   *
   * Functions for removing individual/FORMAT fields
   *
   ******************************************************************************/

  /**
   * @brief Remove an individual field by field name
   *
   * @param tag name of the individual field to remove
   *
   * @note Less efficient than removing using the field index
   */
  VariantBuilder& remove_individual_field(const std::string& tag);

  /**
   * @brief Remove an individual field by field index
   *
   * @param field_index index of the individual field to remove (from a header lookup)
   */
  VariantBuilder& remove_individual_field(const uint32_t field_index);

  /**
   * @brief Remove multiple individual fields by field name
   *
   * @param tags names of the individual fields to remove
   *
   * @note Less efficient than removing using the field indices
   */
  VariantBuilder& remove_individual_fields(const std::vector<std::string>& tags);

  /**
   * @brief Remove multiple individual fields by field index
   *
   * @param field_indices indices of the individual fields to remove (from header lookups)
   */
  VariantBuilder& remove_individual_fields(const std::vector<uint32_t>& field_indices);

  /******************************************************************************
   *
   * Data preparation functions for working with individual fields
   *
   ******************************************************************************/

  /**
   * @brief Get a pre-initialized/padded VariantBuilderMultiSampleVector for use with the more-efficient GT field bulk setters
   *
   * @param num_samples number of samples whose genotypes will be stored in the vector
   * @param max_values_per_sample maximum ploidy (field width)
   *
   * @return a pre-initialized/pre-padded VariantBuilderMultiSampleVector<int32_t> with all samples set to a missing value,
   *         ready for genotypes to be set for non-missing samples
   *
   * @note For usage, see comments to the GT field setters that take a VariantBuilderMultiSampleVector
   */
  VariantBuilderMultiSampleVector<int32_t> get_genotype_multi_sample_vector(const uint32_t num_samples, const uint32_t max_values_per_sample) const;

  /**
   * @brief Get a pre-initialized/padded VariantBuilderMultiSampleVector for use with the more-efficient integer individual field bulk setters
   *
   * @param num_samples number of samples whose values will be stored in the vector
   * @param max_values_per_sample maximum number of values per sample (field width)
   *
   * @return a pre-initialized/pre-padded VariantBuilderMultiSampleVector<int32_t> with all samples set to a missing value,
   *         ready for sample values to be set for non-missing samples
   *
   * @note For usage, see comments to the integer individual field setters that take a VariantBuilderMultiSampleVector
   */
  VariantBuilderMultiSampleVector<int32_t> get_integer_multi_sample_vector(const uint32_t num_samples, const uint32_t max_values_per_sample) const;

  /**
   * @brief Get a pre-initialized/padded VariantBuilderMultiSampleVector for use with the more-efficient float individual field bulk setters
   *
   * @param num_samples number of samples whose values will be stored in the vector
   * @param max_values_per_sample maximum number of values per sample (field width)
   *
   * @return a pre-initialized/pre-padded VariantBuilderMultiSampleVector<float> with all samples set to a missing value,
   *         ready for sample values to be set for non-missing samples
   *
   * @note For usage, see comments to the float individual field setters that take a VariantBuilderMultiSampleVector
   */
  VariantBuilderMultiSampleVector<float> get_float_multi_sample_vector(const uint32_t num_samples, const uint32_t max_values_per_sample) const;


  /******************************************************************************
   *
   * Builder operations: build() and clear()
   *
   ******************************************************************************/

  /**
   * @brief Create a new Variant record using the current state of the builder
   *
   * @return a new Variant object reflecting the field values currently set in this builder
   *
   * @note Can be called multiple times
   * @note It is an error to attempt to build a Variant without having set at least the minimum required
   *       fields (chromosome, alignment start, and reference allele)
   */
  Variant build() const;

  /**
   * @brief Clear all field values in this builder to prepare it for the next build operation
   *
   * @note It's much more efficient to create one VariantBuilder and clear() it between records
   *       instead of creating a new VariantBuilder for each record
   */
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
