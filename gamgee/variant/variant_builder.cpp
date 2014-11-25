#include "variant_builder.h"
#include "missing.h"

#include <algorithm>
#include <stdexcept>

using namespace std;

namespace gamgee {

VariantBuilder::VariantBuilder(const VariantHeader& header) :
  m_header { header.m_header },  // Important: take shared ownership of header rather than make a deep copy
  m_contig {},
  m_start_pos {},
  m_stop_pos {},
  m_qual {},
  m_shared_region { header, true },
  m_individual_region { header, true },
  m_enable_validation { true }
{}

VariantBuilder& VariantBuilder::set_enable_validation(const bool enable_validation) {
  m_enable_validation = enable_validation;
  m_shared_region.set_enable_validation(enable_validation);
  m_individual_region.set_enable_validation(enable_validation);
  return *this;
}

/******************************************************************************
 *
 * Functions for setting core site-level fields
 *
 ******************************************************************************/

VariantBuilder& VariantBuilder::set_chromosome(const uint32_t chromosome) {
  m_contig.set(int32_t(chromosome));
  return *this;
}

VariantBuilder& VariantBuilder::set_chromosome(const std::string& chromosome) {
  // Note: we will validate the contig id (including checking for -1) at build time, if validation is turned on
  m_contig.set(bcf_hdr_id2int(m_header.m_header.get(), BCF_DT_CTG, chromosome.c_str()));
  return *this;
}

VariantBuilder& VariantBuilder::set_alignment_start(const uint32_t alignment_start) {
  m_start_pos.set(alignment_start - 1);
  return *this;
}

VariantBuilder& VariantBuilder::set_alignment_stop(const uint32_t alignment_stop) {
  m_stop_pos.set(alignment_stop - 1);
  return *this;
}

VariantBuilder& VariantBuilder::set_qual(const float qual) {
  m_qual.set(qual);
  return *this;
}

VariantBuilder& VariantBuilder::set_id(const std::string& id) {
  m_shared_region.set_id(id);
  return *this;
}

VariantBuilder& VariantBuilder::set_ref_allele(const std::string& ref_allele) {
  m_shared_region.set_ref_allele(ref_allele);
  return *this;
}

VariantBuilder& VariantBuilder::set_alt_allele(const std::string& alt_allele) {
  m_shared_region.set_alt_allele(alt_allele);
  return *this;
}

VariantBuilder& VariantBuilder::set_alt_alleles(const std::vector<std::string>& alt_alleles) {
  m_shared_region.set_alt_alleles(alt_alleles);
  return *this;
}

VariantBuilder& VariantBuilder::set_filters(const std::vector<std::string>& filters) {
  m_shared_region.set_filters(filters);
  return *this;
}

VariantBuilder& VariantBuilder::set_filters(const std::vector<int32_t>& filters) {
  m_shared_region.set_filters(filters);
  return *this;
}

/******************************************************************************
 *
 * Functions for removing core site-level fields
 *
 ******************************************************************************/

VariantBuilder& VariantBuilder::remove_alignment_stop() {
  m_stop_pos.clear();
  return *this;
}

VariantBuilder& VariantBuilder::remove_qual() {
  m_qual.clear();
  return *this;
}

VariantBuilder& VariantBuilder::remove_id() {
  m_shared_region.remove_id();
  return *this;
}

VariantBuilder& VariantBuilder::remove_alt_alleles() {
  m_shared_region.remove_alt_alleles();
  return *this;
}

VariantBuilder& VariantBuilder::remove_filters() {
  m_shared_region.remove_filters();
  return *this;
}

/******************************************************************************
 *
 * Functions for setting shared/INFO fields
 *
 ******************************************************************************/

VariantBuilder& VariantBuilder::set_integer_shared_field(const std::string& tag, const int32_t value) {
  m_shared_region.set_info_field(tag, value, BCF_HT_INT, 1u);
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_shared_field(const std::string& tag, const std::vector<int32_t>& values) {
  m_shared_region.set_info_field(tag, values, BCF_HT_INT, values.size());
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_shared_field(const uint32_t index, const int32_t value) {
  m_shared_region.set_info_field(index, value, BCF_HT_INT, 1u);
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_shared_field(const uint32_t index, const std::vector<int32_t>& values) {
  m_shared_region.set_info_field(index, values, BCF_HT_INT, values.size());
  return *this;
}

VariantBuilder& VariantBuilder::set_float_shared_field(const std::string& tag, const float value) {
  m_shared_region.set_info_field(tag, value, BCF_HT_REAL, 1u);
  return *this;
}

VariantBuilder& VariantBuilder::set_float_shared_field(const std::string& tag, const std::vector<float>& values) {
  m_shared_region.set_info_field(tag, values, BCF_HT_REAL, values.size());
  return *this;
}

VariantBuilder& VariantBuilder::set_float_shared_field(const uint32_t index, const float value) {
  m_shared_region.set_info_field(index, value, BCF_HT_REAL, 1u);
  return *this;
}

VariantBuilder& VariantBuilder::set_float_shared_field(const uint32_t index, const std::vector<float>& values) {
  m_shared_region.set_info_field(index, values, BCF_HT_REAL, values.size());
  return *this;
}

VariantBuilder& VariantBuilder::set_string_shared_field(const std::string& tag, const std::string& value) {
  m_shared_region.set_info_field(tag, value, BCF_HT_STR, value.length());
  return *this;
}

VariantBuilder& VariantBuilder::set_string_shared_field(const uint32_t index, const std::string& value) {
  m_shared_region.set_info_field(index, value, BCF_HT_STR, value.length());
  return *this;
}

VariantBuilder& VariantBuilder::set_boolean_shared_field(const std::string& tag) {
  m_shared_region.set_info_field(tag, true, BCF_HT_FLAG, 1u);
  return *this;
}

VariantBuilder& VariantBuilder::set_boolean_shared_field(const uint32_t index) {
  m_shared_region.set_info_field(index, true, BCF_HT_FLAG, 1u);
  return *this;
}

/******************************************************************************
 *
 * Functions for removing shared/INFO fields
 *
 ******************************************************************************/

VariantBuilder& VariantBuilder::remove_shared_field(const std::string& tag) {
  m_shared_region.remove_info_field(tag);
  return *this;
}

VariantBuilder& VariantBuilder::remove_shared_field(const uint32_t field_index) {
  m_shared_region.remove_info_field(field_index);
  return *this;
}

VariantBuilder& VariantBuilder::remove_shared_fields(const std::vector<std::string>& tags) {
  for_each(tags.begin(), tags.end(), [this](const std::string& tag){ m_shared_region.remove_info_field(tag); });
  return *this;
}

VariantBuilder& VariantBuilder::remove_shared_fields(const std::vector<uint32_t>& field_indices) {
  for_each(field_indices.begin(), field_indices.end(), [this](const uint32_t field_index){ m_shared_region.remove_info_field(field_index); });
  return *this;
}

/******************************************************************************
 *
 * Functions for setting individual/FORMAT fields in bulk
 * (ie., setting all values for all samples at once)
 *
 ******************************************************************************/

VariantBuilder& VariantBuilder::set_genotypes(const VariantBuilderMultiSampleVector<int32_t>& genotypes_for_all_samples) {
  // Since the user has chosen to pass by lvalue, make a copy before encoding the genotypes
  auto encoded_genotypes = genotypes_for_all_samples;
  Genotype::encode_genotypes(encoded_genotypes);

  // We've made a copy, so we can move the copy into the storage layer
  m_individual_region.bulk_set_genotype_field(m_individual_region.gt_index(), move(encoded_genotypes.get_vector()));
  return *this;
}

VariantBuilder& VariantBuilder::set_genotypes(VariantBuilderMultiSampleVector<int32_t>&& genotypes_for_all_samples) {
  // Encode user's vector directly, since it's been moved in to us
  Genotype::encode_genotypes(genotypes_for_all_samples);

  m_individual_region.bulk_set_genotype_field(m_individual_region.gt_index(), move(genotypes_for_all_samples.get_vector()));
  return *this;
}

VariantBuilder& VariantBuilder::set_genotypes(const std::vector<std::vector<int32_t>>& genotypes_for_all_samples) {
  // Since the user has chosen to pass by lvalue, make a copy before encoding the genotypes
  auto encoded_genotypes = genotypes_for_all_samples;
  Genotype::encode_genotypes(encoded_genotypes);

  // We've made a copy, so we can move the copy into the storage layer
  m_individual_region.bulk_set_genotype_field(m_individual_region.gt_index(), move(encoded_genotypes));
  return *this;
}

VariantBuilder& VariantBuilder::set_genotypes(std::vector<std::vector<int32_t>>&& genotypes_for_all_samples) {
  // Encode user's vector directly, since it's been moved in to us
  Genotype::encode_genotypes(genotypes_for_all_samples);

  m_individual_region.bulk_set_genotype_field(m_individual_region.gt_index(), move(genotypes_for_all_samples));
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_individual_field(const std::string& tag, const VariantBuilderMultiSampleVector<int32_t>& values_for_all_samples) {
  // Ensure that we have an lvalue reference to the values vector so that we make a copy further down the line
  const auto& values_vector = values_for_all_samples.get_vector();
  m_individual_region.bulk_set_integer_field(tag, values_vector);
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_individual_field(const std::string& tag, VariantBuilderMultiSampleVector<int32_t>&& values_for_all_samples) {
  m_individual_region.bulk_set_integer_field(tag, move(values_for_all_samples.get_vector()));
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_individual_field(const std::string& tag, const std::vector<std::vector<int32_t>>& values_for_all_samples) {
  m_individual_region.bulk_set_integer_field(tag, values_for_all_samples);
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_individual_field(const std::string& tag, std::vector<std::vector<int32_t>>&& values_for_all_samples) {
  m_individual_region.bulk_set_integer_field(tag, move(values_for_all_samples));
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_individual_field(const uint32_t field_index, const VariantBuilderMultiSampleVector<int32_t>& values_for_all_samples) {
  // Ensure that we have an lvalue reference to the values vector so that we make a copy further down the line
  const auto& values_vector = values_for_all_samples.get_vector();
  m_individual_region.bulk_set_integer_field(field_index, values_vector);
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_individual_field(const uint32_t field_index, VariantBuilderMultiSampleVector<int32_t>&& values_for_all_samples) {
  m_individual_region.bulk_set_integer_field(field_index, move(values_for_all_samples.get_vector()));
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_individual_field(const uint32_t field_index, const std::vector<std::vector<int32_t>>& values_for_all_samples) {
  m_individual_region.bulk_set_integer_field(field_index, values_for_all_samples);
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_individual_field(const uint32_t field_index, std::vector<std::vector<int32_t>>&& values_for_all_samples) {
  m_individual_region.bulk_set_integer_field(field_index, move(values_for_all_samples));
  return *this;
}

VariantBuilder& VariantBuilder::set_float_individual_field(const std::string& tag, const VariantBuilderMultiSampleVector<float>& values_for_all_samples) {
  // Ensure that we have an lvalue reference to the values vector so that we make a copy further down the line
  const auto& values_vector = values_for_all_samples.get_vector();
  m_individual_region.bulk_set_float_field(tag, values_vector);
  return *this;
}

VariantBuilder& VariantBuilder::set_float_individual_field(const std::string& tag, VariantBuilderMultiSampleVector<float>&& values_for_all_samples) {
  m_individual_region.bulk_set_float_field(tag, move(values_for_all_samples.get_vector()));
  return *this;
}

VariantBuilder& VariantBuilder::set_float_individual_field(const std::string& tag, const std::vector<std::vector<float>>& values_for_all_samples) {
  m_individual_region.bulk_set_float_field(tag, values_for_all_samples);
  return *this;
}

VariantBuilder& VariantBuilder::set_float_individual_field(const std::string& tag, std::vector<std::vector<float>>&& values_for_all_samples) {
  m_individual_region.bulk_set_float_field(tag, move(values_for_all_samples));
  return *this;
}

VariantBuilder& VariantBuilder::set_float_individual_field(const uint32_t field_index, const VariantBuilderMultiSampleVector<float>& values_for_all_samples) {
  // Ensure that we have an lvalue reference to the values vector so that we make a copy further down the line
  const auto& values_vector = values_for_all_samples.get_vector();
  m_individual_region.bulk_set_float_field(field_index, values_vector);
  return *this;
}

VariantBuilder& VariantBuilder::set_float_individual_field(const uint32_t field_index, VariantBuilderMultiSampleVector<float>&& values_for_all_samples) {
  m_individual_region.bulk_set_float_field(field_index, move(values_for_all_samples.get_vector()));
  return *this;
}

VariantBuilder& VariantBuilder::set_float_individual_field(const uint32_t field_index, const std::vector<std::vector<float>>& values_for_all_samples) {
  m_individual_region.bulk_set_float_field(field_index, values_for_all_samples);
  return *this;
}

VariantBuilder& VariantBuilder::set_float_individual_field(const uint32_t field_index, std::vector<std::vector<float>>&& values_for_all_samples) {
  m_individual_region.bulk_set_float_field(field_index, move(values_for_all_samples));
  return *this;
}

VariantBuilder& VariantBuilder::set_string_individual_field(const std::string& tag, const std::vector<std::string>& values_for_all_samples) {
  m_individual_region.bulk_set_string_field(tag, values_for_all_samples);
  return *this;
}

VariantBuilder& VariantBuilder::set_string_individual_field(const std::string& tag, std::vector<std::string>&& values_for_all_samples) {
  m_individual_region.bulk_set_string_field(tag, move(values_for_all_samples));
  return *this;
}

VariantBuilder& VariantBuilder::set_string_individual_field(const uint32_t field_index, const std::vector<std::string>& values_for_all_samples) {
  m_individual_region.bulk_set_string_field(field_index, values_for_all_samples);
  return *this;
}

VariantBuilder& VariantBuilder::set_string_individual_field(const uint32_t field_index, std::vector<std::string>&& values_for_all_samples) {
  m_individual_region.bulk_set_string_field(field_index, move(values_for_all_samples));
  return *this;
}

/******************************************************************************
 *
 * Functions for setting individual/FORMAT fields by sample
 *
 ******************************************************************************/

VariantBuilder& VariantBuilder::set_genotype(const std::string& sample, const std::vector<int32_t>& genotype) {
  // Since the user has passed by lvalue, make a copy before encoding
  auto encoded_genotype = genotype;
  Genotype::encode_genotype(encoded_genotype);

  m_individual_region.set_genotype_field_by_sample(m_individual_region.gt_index(), sample, encoded_genotype.empty() ? nullptr : &(encoded_genotype[0]), encoded_genotype.size());
  return *this;
}

VariantBuilder& VariantBuilder::set_genotype(const std::string& sample, std::vector<int32_t>&& genotype) {
  // Encode user's vector directly, since it's been moved in to us
  Genotype::encode_genotype(genotype);

  m_individual_region.set_genotype_field_by_sample(m_individual_region.gt_index(), sample, genotype.empty() ? nullptr : &(genotype[0]), genotype.size());
  return *this;
}

VariantBuilder& VariantBuilder::set_genotype(const uint32_t sample_index, const std::vector<int32_t>& genotype) {
  // Since the user has passed by lvalue, make a copy before encoding
  auto encoded_genotype = genotype;
  Genotype::encode_genotype(encoded_genotype);

  m_individual_region.set_genotype_field_by_sample(m_individual_region.gt_index(), sample_index, encoded_genotype.empty() ? nullptr : &(encoded_genotype[0]), encoded_genotype.size());
  return *this;
}

VariantBuilder& VariantBuilder::set_genotype(const uint32_t sample_index, std::vector<int32_t>&& genotype) {
  // Encode user's vector directly, since it's been moved in to us
  Genotype::encode_genotype(genotype);

  m_individual_region.set_genotype_field_by_sample(m_individual_region.gt_index(), sample_index, genotype.empty() ? nullptr : &(genotype[0]), genotype.size());
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_individual_field(const std::string& tag, const std::string& sample, const int32_t value) {
  m_individual_region.set_integer_field_by_sample(tag, sample, &value, 1);
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_individual_field(const std::string& tag, const std::string& sample, const std::vector<int32_t>& values) {
  m_individual_region.set_integer_field_by_sample(tag, sample, values.empty() ? nullptr : &(values[0]), values.size());
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_individual_field(const uint32_t field_index, const uint32_t sample_index, const int32_t value) {
  m_individual_region.set_integer_field_by_sample(field_index, sample_index, &value, 1);
  return *this;
}

VariantBuilder& VariantBuilder::set_integer_individual_field(const uint32_t field_index, const uint32_t sample_index, const std::vector<int32_t>& values) {
  m_individual_region.set_integer_field_by_sample(field_index, sample_index, values.empty() ? nullptr : &(values[0]), values.size());
  return *this;
}

VariantBuilder& VariantBuilder::set_float_individual_field(const std::string& tag, const std::string& sample, const float value) {
  m_individual_region.set_float_field_by_sample(tag, sample, &value, 1);
  return *this;
}

VariantBuilder& VariantBuilder::set_float_individual_field(const std::string& tag, const std::string& sample, const std::vector<float>& values) {
  m_individual_region.set_float_field_by_sample(tag, sample, values.empty() ? nullptr : &(values[0]), values.size());
  return *this;
};

VariantBuilder& VariantBuilder::set_float_individual_field(const uint32_t field_index, const uint32_t sample_index, const float value) {
  m_individual_region.set_float_field_by_sample(field_index, sample_index, &value, 1);
  return *this;
}

VariantBuilder& VariantBuilder::set_float_individual_field(const uint32_t field_index, const uint32_t sample_index, const std::vector<float>& values) {
  m_individual_region.set_float_field_by_sample(field_index, sample_index, values.empty() ? nullptr : &(values[0]), values.size());
  return *this;
}

VariantBuilder& VariantBuilder::set_string_individual_field(const std::string& tag, const std::string& sample, const std::string& value) {
  m_individual_region.set_string_field_by_sample(tag, sample, value.empty() ? nullptr : value.c_str(), value.length());
  return *this;
}

VariantBuilder& VariantBuilder::set_string_individual_field(const uint32_t field_index, const uint32_t sample_index, const std::string& value) {
  m_individual_region.set_string_field_by_sample(field_index, sample_index, value.empty() ? nullptr : value.c_str(), value.length());
  return *this;
}

/******************************************************************************
 *
 * Functions for removing individual/FORMAT fields
 *
 ******************************************************************************/

VariantBuilder& VariantBuilder::remove_individual_field(const std::string& tag) {
  m_individual_region.remove_individual_field(tag);
  return *this;
}

VariantBuilder& VariantBuilder::remove_individual_field(const uint32_t field_index) {
  m_individual_region.remove_individual_field(field_index);
  return *this;
}

VariantBuilder& VariantBuilder::remove_individual_fields(const std::vector<std::string>& tags) {
  for_each(tags.begin(), tags.end(), [this](const std::string& tag){ m_individual_region.remove_individual_field(tag); });
  return *this;
}

VariantBuilder& VariantBuilder::remove_individual_fields(const std::vector<uint32_t>& field_indices) {
  for_each(field_indices.begin(), field_indices.end(), [this](const uint32_t field_index){ m_individual_region.remove_individual_field(field_index); });
  return *this;
}


/******************************************************************************
*
* Data preparation functions for working with individual fields
*
******************************************************************************/

VariantBuilderMultiSampleVector<int32_t> VariantBuilder::get_genotype_multi_sample_vector(const uint32_t num_samples, const uint32_t max_values_per_sample) const {
  return VariantBuilderMultiSampleVector<int32_t>{num_samples, max_values_per_sample, -1, bcf_int32_vector_end};
}

VariantBuilderMultiSampleVector<int32_t> VariantBuilder::get_integer_multi_sample_vector(const uint32_t num_samples, const uint32_t max_values_per_sample) const {
  return VariantBuilderMultiSampleVector<int32_t>{num_samples, max_values_per_sample, bcf_int32_missing, bcf_int32_vector_end};
}

VariantBuilderMultiSampleVector<float> VariantBuilder::get_float_multi_sample_vector(const uint32_t num_samples, const uint32_t max_values_per_sample) const {
  auto float_missing = 0.0f; bcf_float_set_missing(float_missing);
  auto float_vector_end = 0.0f; bcf_float_set_vector_end(float_vector_end);
  return VariantBuilderMultiSampleVector<float>{num_samples, max_values_per_sample, float_missing, float_vector_end};
}


/******************************************************************************
 *
 * Builder operations: build() and clear()
 *
 ******************************************************************************/

Variant VariantBuilder::build() const {
  const auto new_variant_body = utils::make_shared_variant(bcf_init1());

  // Always build from scratch for now, until support for building from a starting variant is added
  build_from_scratch(new_variant_body);

  if ( m_enable_validation ) {
    post_build_validation(new_variant_body);
  }

  return Variant{m_header.m_header, new_variant_body};
}

VariantBuilder& VariantBuilder::clear() {
  m_contig.clear();
  m_start_pos.clear();
  m_stop_pos.clear();
  m_qual.clear();

  m_shared_region.clear();
  m_individual_region.clear();

  return *this;
}

/******************************************************************************
 *
 * Private functionality
 *
 ******************************************************************************/

void VariantBuilder::build_from_scratch(const std::shared_ptr<bcf1_t>& new_variant_body) const {
  // Missing rid/pos are errors that will be caught in post-build validation (if validation is turned on)
  new_variant_body->rid = m_contig.is_set() ? m_contig.field_value() : missing_values::int32;
  new_variant_body->pos = m_start_pos.is_set() ? m_start_pos.field_value() : missing_values::int32;

  if ( m_qual.is_set() ) {
    new_variant_body->qual = m_qual.field_value();
  }
  else {
    bcf_float_set_missing(new_variant_body->qual);
  }

  // Set rlen to the reference block size if alignment stop was set, otherwise set it to the ref allele length.
  // If m_start_pos is not set, we'll get an error in the post-validation stage.
  new_variant_body->rlen = m_stop_pos.is_set() ? (m_stop_pos.field_value() - m_start_pos.field_value() + 1) :
                                                 int32_t(m_shared_region.ref_allele_length());
  new_variant_body->n_allele = 1 + m_shared_region.num_alt_alleles();

  // Shared region (always encoded, since at a minimum the ref allele will be present)
  new_variant_body->n_info = m_shared_region.num_present_info_fields();
  auto shared_buffer = utils::initialize_htslib_buffer(m_shared_region.estimate_total_size());
  m_shared_region.encode_into(&shared_buffer);
  new_variant_body->shared = shared_buffer;

  // Individual region (conditionally encoded)
  new_variant_body->n_sample = m_header.n_samples();
  new_variant_body->n_fmt = m_individual_region.num_present_fields();

  if ( m_individual_region.num_present_fields() > 0 ) {
    auto indiv_buffer = utils::initialize_htslib_buffer(m_individual_region.estimate_total_size());
    m_individual_region.encode_into(&indiv_buffer);
    new_variant_body->indiv = indiv_buffer;
  }
  else {
    new_variant_body->indiv = {0, 0, 0};
  }
}

void VariantBuilder::post_build_validation(const std::shared_ptr<bcf1_t>& new_variant_body) const {
  // Note that we only do validation of core (non-data-region) fields here.
  // The data region fields are validated in the VariantBuilderSharedRegion and
  // VariantBuilderIndividualRegion/VariantBuilderIndividualField classes.

  // Chromosome field must be present
  if ( new_variant_body->rid == missing_values::int32 ) {
    throw logic_error("Missing required chromosome field");
  }

  // Chromosome must be in header sequence dictionary
  if ( new_variant_body->rid < 0 || new_variant_body->rid >= m_header.m_header->n[BCF_DT_CTG] || m_header.m_header->id[BCF_DT_CTG][new_variant_body->rid].key == nullptr ) {
    throw invalid_argument(string{"Chromosome with index "} + to_string(new_variant_body->rid) + " not found in header sequence dictionary");
  }

  // Alignment start field must be present
  if ( new_variant_body->pos == missing_values::int32 ) {
    throw logic_error("Missing required alignment start field");
  }

  // rlen will be negative or 0 at this point if alignment stop was < alignment start
  if ( new_variant_body->rlen <= 0 ) {
    throw logic_error("Alignment stop must be >= alignment start");
  }
}

}
