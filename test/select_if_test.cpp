#include <boost/test/unit_test.hpp>
#include <boost/dynamic_bitset.hpp>
#include "variant/variant_reader.h"
#include "variant/variant.h"


using namespace std;
using namespace gamgee;
using namespace boost;

constexpr auto GQ_THRESH = 35;

BOOST_AUTO_TEST_CASE( select_if_individual_fields ) {
  const auto actual_gq_selects = vector<dynamic_bitset<>>{
    dynamic_bitset<>(string("111")),  //Beware that in setting dynamic_bitset that way
    dynamic_bitset<>(string("110")),  //the bit order is inversed:  Here it is actually 011
    dynamic_bitset<>(string("100")),
    dynamic_bitset<>(string("010")),
    dynamic_bitset<>(string("011"))};
  const auto actual_pl_selects = vector<dynamic_bitset<>>{
    dynamic_bitset<>(string("010")),
    dynamic_bitset<>(string("000")),
    dynamic_bitset<>(string("001")),
    dynamic_bitset<>(string("111")),
    dynamic_bitset<>(string("100"))};
  for (const auto& filename : {"testdata/test_variants_02.vcf"}) {
    const auto reader = SingleVariantReader{filename};
    auto record_idx = 0u;
    for (const auto& record : SingleVariantReader{filename}) {
      const auto g_quals = record.integer_individual_field("GQ");
      const auto p_likes = record.integer_individual_field("PL");
      const auto high_gq = [](const IndividualFieldValue<int32_t>& x) { return x[0] >= GQ_THRESH; };
      const auto hom_ref = [](const IndividualFieldValue<int32_t>& x) { return x[0] == 0; };
      const auto comput_gq_select = Variant::select_if(g_quals.begin(), g_quals.end(), high_gq);
      const auto comput_pl_select = Variant::select_if(p_likes.begin(), p_likes.end(), hom_ref);
      BOOST_CHECK_EQUAL(comput_gq_select.size(), 3u);
      BOOST_CHECK_EQUAL(comput_pl_select.size(), 3u);
      const auto actual_gq_select = actual_gq_selects[record_idx];
      const auto actual_pl_select = actual_pl_selects[record_idx];
      BOOST_CHECK(comput_pl_select == actual_pl_select);
      ++record_idx;
    }
    BOOST_CHECK_EQUAL(record_idx, 5u);
  }
}

BOOST_AUTO_TEST_CASE( select_if_shared_field ) {
  const auto truth_an_counts = std::vector<uint32_t>{1,1,1,1,1,1,1};
  const auto truth_af_counts = std::vector<uint32_t>{0,0,0,0,1,0,1};
  auto truth_index = 0u;
  for (const auto& record : SingleVariantReader{"testdata/test_variants.vcf"}) {
    const auto an = record.integer_shared_field("AN");
    const auto r1 = Variant::select_if(an.begin(), an.end(), [](const int v) { return v == 6; });
    BOOST_CHECK_EQUAL(r1.count(), truth_an_counts[truth_index]);
    const auto af = record.float_shared_field("AF");
    const auto r2 = Variant::select_if(af.begin(), af.end(), [](const float v) { return v < 0.5; });
    BOOST_CHECK_EQUAL(r2.count(), truth_af_counts[truth_index++]);
  }
}
