#include <boost/test/unit_test.hpp>
#include <boost/dynamic_bitset.hpp>
#include "variant_reader.h"
#include "variant.h"


using namespace std;
using namespace gamgee;
using namespace boost;

constexpr auto GQ_THRESH = 35;

BOOST_AUTO_TEST_CASE( select_if ) {
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
      const auto g_quals = record.genotype_quals();
      const auto p_likes = record.phred_likelihoods();
      const auto high_gq = [](const IndividualFieldValue<int32_t>& x) { return x[0] >= GQ_THRESH; };
      const auto hom_ref = [](const IndividualFieldValue<int32_t>& x) { return x[0] == 0; };
      const auto comput_gq_select = Variant::select_if(g_quals.begin(), g_quals.end(), high_gq);
      const auto comput_pl_select = Variant::select_if(p_likes.begin(), p_likes.end(), hom_ref);
      BOOST_CHECK_EQUAL(comput_gq_select.size(), 3);
      BOOST_CHECK_EQUAL(comput_pl_select.size(), 3);
      const auto actual_gq_select = actual_gq_selects[record_idx];
      const auto actual_pl_select = actual_pl_selects[record_idx];
      BOOST_CHECK(comput_pl_select == actual_pl_select);
      ++record_idx;
    }
    BOOST_CHECK_EQUAL(record_idx, 5u);
  }
}


