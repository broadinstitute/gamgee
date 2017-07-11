#include "test_utils.h"
#include "variant/variant_header_builder.h"
#include "missing.h"

#include <set>

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

template<class T>
void check_fields(T actual, T truth) {
  std::sort(actual.begin(), actual.end());
  std::sort(truth.begin(), truth.end());
  BOOST_CHECK_EQUAL(actual.size(), truth.size());
  BOOST_CHECK(actual == truth);
}

const auto samples     = vector<string>{"S1", "S292", "S30034"};
const auto sample_indices = vector<int32_t>{0, 1, 2};
const auto chromosomes = vector<string>{"chr1", "chr2", "chr3", "chr4"};
const auto filters     = vector<string>{"LOW_QUAL", "PASS", "VQSR_FAILED"};
const auto filter_indices = vector<int32_t>{1, 0, 2};
const auto shareds       = vector<string>{"DP", "MQ", "RankSum"};
const auto shareds_indices = vector<int32_t>{3,4,5};  // looks arbitrary but these are the indices of the shared fields because the filters get 0, 1 and 2.
const auto individuals     = vector<string>{"GQ", "PL", "DP"};
const auto individuals_indices = vector<int32_t>{6,7,3}; // the last index gets the same number as the info index. Weird, but that's how htslib deals with this.

void variant_header_builder_checks(const VariantHeader& vh) {
  check_fields(vh.chromosomes(), chromosomes);
  check_fields(vh.samples(), samples);
  check_fields(vh.filters(), filters);
  check_fields(vh.shared_fields(), shareds);
  check_fields(vh.individual_fields(), individuals);

  BOOST_CHECK(vh.has_filter("PASS"));
  BOOST_CHECK(vh.has_filter(vh.field_index("PASS")));
  BOOST_CHECK(vh.has_filter("LOW_QUAL"));
  BOOST_CHECK(vh.has_filter(vh.field_index("LOW_QUAL")));
  BOOST_CHECK(vh.has_filter("VQSR_FAILED"));
  BOOST_CHECK(vh.has_filter(vh.field_index("VQSR_FAILED")));
  BOOST_CHECK(! vh.has_filter("BLAH"));
  BOOST_CHECK(! vh.has_filter(vh.field_index("BLAH")));
  BOOST_CHECK(! vh.has_filter("DP"));
  BOOST_CHECK(! vh.has_filter(vh.field_index("DP")));
  BOOST_CHECK(! vh.has_filter("GQ"));
  BOOST_CHECK(! vh.has_filter(vh.field_index("GQ")));
  BOOST_CHECK(! vh.has_filter(20));

  BOOST_CHECK(vh.has_shared_field("DP"));
  BOOST_CHECK(vh.has_shared_field(vh.field_index("DP")));
  BOOST_CHECK(vh.has_shared_field("MQ"));
  BOOST_CHECK(vh.has_shared_field(vh.field_index("MQ")));
  BOOST_CHECK(vh.has_shared_field("RankSum"));
  BOOST_CHECK(vh.has_shared_field(vh.field_index("RankSum")));
  BOOST_CHECK(! vh.has_shared_field("BLAH"));
  BOOST_CHECK(! vh.has_shared_field(vh.field_index("BLAH")));
  BOOST_CHECK(! vh.has_shared_field("LOW_QUAL"));
  BOOST_CHECK(! vh.has_shared_field(vh.field_index("LOW_QUAL")));
  BOOST_CHECK(! vh.has_shared_field("GQ"));
  BOOST_CHECK(! vh.has_shared_field(vh.field_index("GQ")));
  BOOST_CHECK(! vh.has_shared_field(20));

  BOOST_CHECK(vh.has_individual_field("GQ"));
  BOOST_CHECK(vh.has_individual_field(vh.field_index("GQ")));
  BOOST_CHECK(vh.has_individual_field("PL"));
  BOOST_CHECK(vh.has_individual_field(vh.field_index("PL")));
  BOOST_CHECK(vh.has_individual_field("DP"));
  BOOST_CHECK(vh.has_individual_field(vh.field_index("DP")));
  BOOST_CHECK(! vh.has_individual_field("BLAH"));
  BOOST_CHECK(! vh.has_individual_field(vh.field_index("BLAH")));
  BOOST_CHECK(! vh.has_individual_field("MQ"));
  BOOST_CHECK(! vh.has_individual_field(vh.field_index("MQ")));
  BOOST_CHECK(! vh.has_individual_field("LOW_QUAL"));
  BOOST_CHECK(! vh.has_individual_field(vh.field_index("LOW_QUAL")));
  BOOST_CHECK(! vh.has_individual_field(20));

  // Also test API functions that take an explicit field category as a parameter (BCF_HL_INFO, etc.)
  BOOST_CHECK(vh.has_field("MQ", BCF_HL_INFO));
  BOOST_CHECK(vh.has_field(vh.field_index("MQ"), BCF_HL_INFO));
  BOOST_CHECK(! vh.has_field("BLAH", BCF_HL_INFO));
  BOOST_CHECK(! vh.has_field(vh.field_index("BLAH"), BCF_HL_INFO));
  BOOST_CHECK(vh.has_field("GQ", BCF_HL_FMT));
  BOOST_CHECK(vh.has_field(vh.field_index("GQ"), BCF_HL_FMT));
  BOOST_CHECK(! vh.has_field("BLAH", BCF_HL_FMT));
  BOOST_CHECK(! vh.has_field(vh.field_index("BLAH"), BCF_HL_FMT));

  // Test type-querying functions (note: these assume that you've already checked field existence)
  BOOST_CHECK_EQUAL(vh.shared_field_type("MQ"), BCF_HT_INT);
  BOOST_CHECK_EQUAL(vh.shared_field_type(vh.field_index("MQ")), BCF_HT_INT);
  BOOST_CHECK_EQUAL(vh.field_type("MQ", BCF_HL_INFO), BCF_HT_INT);
  BOOST_CHECK_EQUAL(vh.field_type(vh.field_index("MQ"), BCF_HL_INFO), BCF_HT_INT);
  BOOST_CHECK_EQUAL(vh.individual_field_type("PL"), BCF_HT_REAL);
  BOOST_CHECK_EQUAL(vh.individual_field_type(vh.field_index("PL")), BCF_HT_REAL);
  BOOST_CHECK_EQUAL(vh.field_type("PL", BCF_HL_FMT), BCF_HT_REAL);
  BOOST_CHECK_EQUAL(vh.field_type(vh.field_index("PL"), BCF_HL_FMT), BCF_HT_REAL);

  BOOST_CHECK(vh.has_sample("S1"));
  BOOST_CHECK(vh.has_sample(vh.sample_index("S1")));
  BOOST_CHECK(vh.has_sample("S292"));
  BOOST_CHECK(vh.has_sample(vh.sample_index("S292")));
  BOOST_CHECK(vh.has_sample("S30034"));
  BOOST_CHECK(vh.has_sample(vh.sample_index("S30034")));
  BOOST_CHECK(! vh.has_sample("BLAH"));
  BOOST_CHECK(! vh.has_sample(vh.sample_index("BLAH")));
  BOOST_CHECK(! vh.has_sample("DP"));
  BOOST_CHECK(! vh.has_sample(vh.sample_index("DP")));
  BOOST_CHECK(! vh.has_sample(20));

  for (auto i = 0u; i != filters.size(); ++i)
    BOOST_CHECK_EQUAL(filter_indices[i], vh.field_index(filters[i]));
  for (auto i = 0u; i != shareds.size(); ++i)
    BOOST_CHECK_EQUAL(shareds_indices[i], vh.field_index(shareds[i]));
  for (auto i = 0u; i != individuals.size(); ++i)
    BOOST_CHECK_EQUAL(individuals_indices[i], vh.field_index(individuals[i]));
  for (auto i = 0u; i != samples.size(); ++i)
    BOOST_CHECK_EQUAL(sample_indices[i], vh.sample_index(samples[i]));

  BOOST_CHECK(missing(vh.field_index("MISSING")));
  BOOST_CHECK(missing(vh.sample_index("MISSING")));
}

VariantHeaderBuilder simple_builder() {
  auto builder = VariantHeaderBuilder{};
  builder.add_source("Gamgee api test");
  builder.advanced_add_arbitrary_line("##unused=<XX=AA,Description=\"Unused generic\">");
  for (const auto& chromosome : chromosomes)
    builder.add_chromosome(chromosome, "234");
  for (const auto& sample : samples)
    builder.add_sample(sample);
  for (const auto& filter : filters)
    builder.add_filter(filter, "anything", "deer=vanison");
  for (const auto& shared : shareds)
    builder.add_shared_field(shared, "43", "Integer", "something", "the_source", "3.4", "cow=beef");
  for (const auto& individual : individuals)
    builder.add_individual_field(individual, "13", "Float", "nothing", "goat=shank");

  return builder;
}

BOOST_AUTO_TEST_CASE( variant_header_builder_simple_building ) {
  auto builder = simple_builder();
  variant_header_builder_checks(builder.build());
}

BOOST_AUTO_TEST_CASE( variant_header_builder_reuse ) {
  auto builder = simple_builder();

  const auto header1 = builder.build();
  variant_header_builder_checks(header1);

  builder.add_filter("BLAH", "anything", "deer=vanison");
  builder.add_shared_field("BLAH", "43", "Integer", "something", "the_source", "3.4", "cow=beef");
  builder.add_individual_field("BLAH", "13", "Float", "nothing", "goat=shank");
  const auto header2 = builder.build();

  BOOST_CHECK(header2.has_filter("BLAH"));
  BOOST_CHECK(header2.has_shared_field("BLAH"));
  BOOST_CHECK(header2.has_individual_field("BLAH"));

  // these check for the absence of the "BLAH" fields added above
  variant_header_builder_checks(header1);
}

BOOST_AUTO_TEST_CASE( variant_header_builder_one_time_build ) {
  auto builder = simple_builder();

  const auto header = builder.one_time_build();
  variant_header_builder_checks(header);
}

BOOST_AUTO_TEST_CASE( variant_header_builder_move ) {
  auto builder1 = simple_builder();
  auto builder2 = check_move_constructor(builder1);
  variant_header_builder_checks(builder2.build());
}

BOOST_AUTO_TEST_CASE( variant_header_builder_chained ) {
  auto builder = VariantHeaderBuilder{};
  builder
    .add_source("Gamgee api test")
    .advanced_add_arbitrary_line("##unused=<XX=AA,Description=\"Unused generic\">")
    .add_chromosome("chr1", "234")
    .add_chromosome("chr2", "234")
    .add_chromosome("chr3", "234");

  // verify that the chain can be broken and resumed
  builder
    .add_chromosome("chr4", "234")
    .add_sample("S1")
    .add_sample("S292")
    .add_sample("S30034")
    .add_filter("LOW_QUAL", "anything", "deer=vanison")
    .add_filter("PASS", "anything", "deer=vanison")
    .add_filter("VQSR_FAILED", "anything", "deer=vanison")
    .add_shared_field("DP", "43", "Integer", "something", "the_source", "3.4", "cow=beef")
    .add_shared_field("MQ", "43", "Integer", "something", "the_source", "3.4", "cow=beef")
    .add_shared_field("RankSum", "43", "Integer", "something", "the_source", "3.4", "cow=beef")
    .add_individual_field("GQ", "13", "Float", "nothing", "goat=shank")
    .add_individual_field("PL", "13", "Float", "nothing", "goat=shank")
    .add_individual_field("DP", "13", "Float", "nothing", "goat=shank");

  for (const auto& filter : filters)
    builder.add_filter(filter, "anything", "deer=vanison");
  for (const auto& shared : shareds)
    builder.add_shared_field(shared, "43", "Integer", "something", "the_source", "3.4", "cow=beef");
  for (const auto& individual : individuals)
    builder.add_individual_field(individual, "13", "Float", "nothing", "goat=shank");

  variant_header_builder_checks(builder.build());
}

const auto merge_file_1 = "testdata/var_hdr_merge/test1.vcf";
const auto merge_file_2 = "testdata/var_hdr_merge/test2.vcf";
const auto merge_file_3 = "testdata/var_hdr_merge/test3.vcf";

BOOST_AUTO_TEST_CASE( variant_header_merge_test ) {
  auto header1 = SingleVariantReader{merge_file_1}.header();
  auto header2 = SingleVariantReader{merge_file_2}.header();
  auto header3 = SingleVariantReader{merge_file_3}.header();
  BOOST_CHECK(header1 != header2);
  BOOST_CHECK(header1 != header3);
  BOOST_CHECK(header2 != header3);

  // header 3 = the contents of header 1 and header 2 combined

  auto builder = VariantHeaderBuilder{header1};
  builder.merge(header2);

  auto built = builder.build();
  BOOST_CHECK(built != header1);
  BOOST_CHECK(built != header2);
  BOOST_CHECK(built == header3);
}

BOOST_AUTO_TEST_CASE( variant_header_file_and_construct ) {
  auto header1 = SingleVariantReader{merge_file_1}.header();
  auto header2 = SingleVariantReader{merge_file_2}.header();
  auto header3 = SingleVariantReader{merge_file_3}.header();

  // start with header1 contents
  auto builder_from_header = VariantHeaderBuilder{header1};

  // add header2 contents
  builder_from_header
    .add_chromosome("20", "64000000")
    .add_chromosome("22", "120000000")
    .add_sample("SAMPLE2")
    .add_sample("SAMPLE3")
    .add_filter("LOW_QUAL", "Low quality call", "")
    .add_filter("MISSED", "Missed by the variant caller", "")
    .add_shared_field("AN", "1", "Integer", "Total number of alleles in called genotypes", "", "", "")
    .add_shared_field("VALIDATED", "0", "Flag", "Validated By Follow-up Experiment", "", "", "")
    .add_individual_field("GQ", "1", "Integer", "Genotype quality", "")
    .add_individual_field("PL", "G", "Integer", "Phred scaled relative Likelihoods of the genotypes", "");

  BOOST_CHECK(builder_from_header.build() == header3);

  // start with header1 contents
  auto builder_from_scratch = VariantHeaderBuilder{};
  builder_from_scratch
    .add_chromosome("1", "300000000")
    .add_sample("SAMPLE1")
    .add_filter("PASS", "All filters passed", "")
    // need to escape the quotes because the text contains commas
    .add_shared_field("AF", "A", "Float", "\"Allele Frequency, for each ALT allele, in the same order as listed\"", "", "", "")
    .add_individual_field("GT", "1", "String", "Genotype", "");

  // add header2 contents
  builder_from_scratch.merge(header2);

  BOOST_CHECK(builder_from_scratch.build() == header3);
}

const auto filter_truth = vector<bool> { true, true, true, false, false, false, false, false };
const auto shared_truth = vector<bool> { false, false, false, true, true, true, false, false };
const auto individual_truth = vector<bool> { false, false, false, true, false, false, true, true };

// remove duplicate names: assign the same index

void add_to_index_names(const vector<string>& new_names, set<string>& index_name_set, vector<string>& index_name_vector) {
  for (const auto& name : new_names) {
    const auto& finder = index_name_set.find(name);
    if (finder == index_name_set.end()) {
      index_name_set.insert(name);
      index_name_vector.push_back(name);
    }
  }
}

BOOST_AUTO_TEST_CASE( variant_header_field_index_iteration ) {
  auto header = simple_builder().build();

  auto index_name_set = set<string>{};
  auto index_name_vector = vector<string>{};

  // the PASS filter is always index 0
  add_to_index_names({"PASS"}, index_name_set, index_name_vector);
  add_to_index_names(filters, index_name_set, index_name_vector);
  add_to_index_names(shareds, index_name_set, index_name_vector);
  add_to_index_names(individuals, index_name_set, index_name_vector);

  const auto& name_truth = index_name_vector;
  const auto end_truth = index_name_vector.size();

  BOOST_CHECK_EQUAL(header.field_index_end(), end_truth);

  for (auto idx = 0u; idx < header.field_index_end(); ++idx) {
    BOOST_CHECK_EQUAL(header.has_filter(idx), filter_truth[idx]);
    BOOST_CHECK_EQUAL(header.has_shared_field(idx), shared_truth[idx]);
    BOOST_CHECK_EQUAL(header.has_individual_field(idx), individual_truth[idx]);
    BOOST_CHECK_EQUAL(header.get_field_name(idx), name_truth[idx]);
  }
}
