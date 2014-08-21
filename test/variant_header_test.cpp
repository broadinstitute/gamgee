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

BOOST_AUTO_TEST_CASE( variant_header_builder_simple_building ) {
  const auto samples     = vector<string>{"S1", "S292", "S30034"};
  const auto chromosomes = vector<string>{"chr1", "chr2", "chr3", "chr4"};
  const auto filters     = vector<string>{"LOW_QUAL", "PASS", "VQSR_FAILED"};
  const auto shareds       = vector<string>{"DP", "MQ", "RankSum"};
  const auto shareds_indices = vector<uint32_t>{3,4,5};  // looks arbitrary but these are the indices of the shared fields because the filters get 0, 1 and 2.
  const auto individuals     = vector<string>{"GQ", "PL", "DP"};
  const auto individuals_indices = vector<uint32_t>{6,7,3}; // the last index gets the same number as the info index. Weird, but that's how htslib deals with this.
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
  const auto vh = builder.build();
  check_fields(vh.chromosomes(), chromosomes);
  check_fields(vh.samples(), samples);
  check_fields(vh.filters(), filters);
  check_fields(vh.shared_fields(), shareds);
  check_fields(vh.individual_fields(), individuals);
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

BOOST_AUTO_TEST_CASE( variant_header_move_and_copy_constructor ) {
  auto builder = VariantHeaderBuilder{};
  builder.add_sample("S1");
  auto h0 = builder.build();
  auto copies = check_copy_constructor(h0);
  auto c2 = get<2>(copies);
  BOOST_CHECK_EQUAL_COLLECTIONS(h0.samples().begin(), h0.samples().end(), get<0>(copies).samples().begin(), get<0>(copies).samples().end());
  BOOST_CHECK_EQUAL_COLLECTIONS(h0.samples().begin(), h0.samples().end(), get<1>(copies).samples().begin(), get<1>(copies).samples().end());
  BOOST_CHECK_EQUAL_COLLECTIONS(h0.samples().begin(), h0.samples().end(), c2.samples().begin(), c2.samples().end());
  auto m1 = check_move_constructor(get<1>(copies));
  BOOST_CHECK_EQUAL_COLLECTIONS(h0.samples().begin(), h0.samples().end(), m1.samples().begin(), m1.samples().end());
  // can't modify a variant header... so this is it for the test.
}

