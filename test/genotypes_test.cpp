#include <boost/test/unit_test.hpp>
#include <boost/dynamic_bitset.hpp>
#include "variant_reader.h"
#include "variant.h"
#include "genotype.h"


using namespace std;
using namespace gamgee;
using namespace boost;


void select_if_test(const vector<dynamic_bitset<>>& truth, const std::string& filename, const std::function<bool (const gamgee::Genotype& )> pred){
  auto record_idx = 0u;
  for (const auto& record : SingleVariantReader{filename}) {
    const auto genotype = record.genotypes();
    const auto result = Variant::select_if(genotype.begin(), genotype.end(), pred );

    BOOST_CHECK(result == truth[record_idx]);
    ++record_idx;
  }
}

const auto diploid= "testdata/test_variants.vcf";
const auto multi_ploidy = "testdata/test_variants_alternate_ploidy.vcf";

bool hom_ref(const gamgee::Genotype& g) { return g.hom_ref(); };
bool hom_var(const gamgee::Genotype& g) { return g.hom_var(); };

BOOST_AUTO_TEST_CASE( hom_ref_test ){
  const auto truth = vector<dynamic_bitset<>>{
    dynamic_bitset<>(string("010")),
    dynamic_bitset<>(string("010")),
    dynamic_bitset<>(string("010")),
    dynamic_bitset<>(string("010")),
    dynamic_bitset<>(string("010"))};

  select_if_test(truth, diploid, hom_ref  );
}

BOOST_AUTO_TEST_CASE( hom_ref_alternate_ploidy_test ){
  const auto truth = vector<dynamic_bitset<>>{
    dynamic_bitset<>(string("010")),
    dynamic_bitset<>(string("010")),
    dynamic_bitset<>(string("010"))};

  select_if_test(truth, multi_ploidy , hom_ref );
}

BOOST_AUTO_TEST_CASE( hom_var_test ){
  const auto truth = vector<dynamic_bitset<>>{
    dynamic_bitset<>(string("100")),
    dynamic_bitset<>(string("100")),
    dynamic_bitset<>(string("100")),
    dynamic_bitset<>(string("100")),
    dynamic_bitset<>(string("100"))};

  select_if_test(truth, diploid, hom_var);
}

BOOST_AUTO_TEST_CASE( hom_var_alternate_ploidy_test ){
  const auto truth = vector<dynamic_bitset<>>{
    dynamic_bitset<>(string("101")),
    dynamic_bitset<>(string("100")),
    dynamic_bitset<>(string("100"))};

  select_if_test(truth, multi_ploidy, hom_var);
}

BOOST_AUTO_TEST_CASE( het_alternate_ploidy_test ){
  const auto truth = vector<dynamic_bitset<>>{
    dynamic_bitset<>(string("000")),
    dynamic_bitset<>(string("001")),
    dynamic_bitset<>(string("000"))};

  select_if_test(truth, multi_ploidy, [](const auto& g) {return g.het(); });
}


