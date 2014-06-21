#include <boost/test/unit_test.hpp>
#include "variant_header_builder.h"

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
  const auto samples = vector<string>{"S1", "S292", "S30034"}; 
  const auto contigs = vector<string>{"chr1", "chr2", "chr3", "chr4"}; 
  const auto filters = vector<string>{"LOW_QUAL", "PASS", "VQSR_FAILED"};
  const auto infos   = vector<string>{"DP", "MQ", "RankSum"};
  const auto formats = vector<string>{"GQ", "PL", "DP"};
  auto builder = VariantHeaderBuilder{};
  builder.add_source("Gamgee api test");
  builder.advanced_add_arbitrary_line("##unused=<XX=AA,Description=\"Unused generic\">");
  for (const auto& contig : contigs) 
    builder.add_contig(contig, "234");
  for (const auto& sample : samples)
    builder.add_sample(sample);
  for (const auto& filter : filters) 
    builder.add_filter(filter, "anything", "deer=vanison");
  for (const auto& info : infos) 
    builder.add_info_field(info, "43", "Integer", "something", "the_source", "3.4", "cow=beef");
  for (const auto& format : formats) 
    builder.add_format_field(format, "13", "Float", "nothing", "goat=shank");
  const auto vh = builder.build();
  check_fields(vh.contigs(), contigs);
  check_fields(vh.samples(), samples);
  check_fields(vh.filters(), filters);
  check_fields(vh.info_fields(), infos);
  check_fields(vh.format_fields(), formats);
  BOOST_CHECK(vh.has_format_field("GQ") == true);
  BOOST_CHECK(vh.has_format_field("DP") == true);
  BOOST_CHECK(vh.has_format_field("BLAH") == false);
}

