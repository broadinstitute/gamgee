#include <boost/test/unit_test.hpp>
#include "test_utils.h"
#include "variant_header_builder.h"
#include "missing.h"

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
const auto chromosomes = vector<string>{"chr1", "chr2", "chr3", "chr4"};
const auto filters     = vector<string>{"LOW_QUAL", "PASS", "VQSR_FAILED"};
const auto shareds       = vector<string>{"DP", "MQ", "RankSum"};
const auto shareds_indices = vector<uint32_t>{3,4,5};  // looks arbitrary but these are the indices of the shared fields because the filters get 0, 1 and 2.
const auto individuals     = vector<string>{"GQ", "PL", "DP"};
const auto individuals_indices = vector<uint32_t>{6,7,3}; // the last index gets the same number as the info index. Weird, but that's how htslib deals with this.

void variant_header_builder_checks(const VariantHeader& vh) {
  check_fields(vh.chromosomes(), chromosomes);
  check_fields(vh.samples(), samples);
  check_fields(vh.filters(), filters);
  check_fields(vh.shared_fields(), shareds);
  check_fields(vh.individual_fields(), individuals);
  BOOST_CHECK(vh.has_filter("LOW_QUAL") == true);
  BOOST_CHECK(vh.has_filter("VQSR_FAILED") == true);
  BOOST_CHECK(vh.has_filter("BLAH") == false);
  BOOST_CHECK(vh.has_shared_field("DP") == true);
  BOOST_CHECK(vh.has_shared_field("MQ") == true);
  BOOST_CHECK(vh.has_shared_field("BLAH") == false);
  BOOST_CHECK(vh.has_individual_field("GQ") == true);
  BOOST_CHECK(vh.has_individual_field("DP") == true);
  BOOST_CHECK(vh.has_individual_field("BLAH") == false);
  for (auto i = 0u; i != shareds.size(); ++i)
    BOOST_CHECK_EQUAL(shareds_indices[i], vh.field_index(shareds[i]));
  for (auto i = 0u; i != individuals.size(); ++i)
    BOOST_CHECK_EQUAL(individuals_indices[i], vh.field_index(individuals[i]));
  BOOST_CHECK(missing(vh.field_index("MISSING")));
  BOOST_CHECK(missing(vh.field_index("MISSING")));
  BOOST_CHECK_EQUAL(vh.field_index("PASS"), 0);
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

  for (const auto& sample : samples)
    builder.add_sample(sample);
  for (const auto& filter : filters)
    builder.add_filter(filter, "anything", "deer=vanison");
  for (const auto& shared : shareds)
    builder.add_shared_field(shared, "43", "Integer", "something", "the_source", "3.4", "cow=beef");
  for (const auto& individual : individuals)
    builder.add_individual_field(individual, "13", "Float", "nothing", "goat=shank");

  variant_header_builder_checks(builder.build());
}
