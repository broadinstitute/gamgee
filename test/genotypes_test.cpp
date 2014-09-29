#include <boost/test/unit_test.hpp>

#include "variant_reader.h"
#include "variant.h"
#include "genotype.h"

#include <boost/dynamic_bitset.hpp>

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

BOOST_AUTO_TEST_CASE( het_snp_is_true ) {
  const auto rec = *(SingleVariantReader{"testdata/test_variants.bcf"}.begin());
  const auto mask = rec.allele_mask();
  BOOST_CHECK(rec.genotypes()[0].snp(mask));
  BOOST_CHECK(!rec.genotypes()[0].insertion(mask));
  BOOST_CHECK(!rec.genotypes()[0].deletion(mask));
}

BOOST_AUTO_TEST_CASE( het_insertion_is_true ) {
  auto it = SingleVariantReader{"testdata/test_variants.bcf"}.begin();
  ++it; ++it; ++it;
  const auto rec = *it;
  const auto mask = rec.allele_mask();
  BOOST_CHECK(rec.genotypes()[0].insertion(mask));
  BOOST_CHECK(!rec.genotypes()[0].snp(mask));
  BOOST_CHECK(!rec.genotypes()[0].deletion(mask));
}

BOOST_AUTO_TEST_CASE( het_deletion_is_true ) {
  auto it = SingleVariantReader{"testdata/test_variants.bcf"}.begin();
  ++it; ++it; 
  const auto rec = *it;
  const auto mask = rec.allele_mask();
  BOOST_CHECK(rec.genotypes()[0].deletion(mask));
  BOOST_CHECK(!rec.genotypes()[0].insertion(mask));
  BOOST_CHECK(!rec.genotypes()[0].snp(mask));
}

BOOST_AUTO_TEST_CASE( het_indel_is_true ) {
  const auto test_rec = [](const Variant& r, const bool t){
    auto mask = r.allele_mask();
    BOOST_CHECK_EQUAL(r.genotypes()[0].indel(mask), t);
  };
  const auto truth = vector<bool>{false, false, true, true, true, true};
  auto i = 0;
  for (const auto& rec : SingleVariantReader{"testdata/test_variants.bcf"}) {
    test_rec(rec, truth[i]);
    ++i;
  }
}

BOOST_AUTO_TEST_CASE( complex_test ) {
  auto test_rec = [](const Variant& r, const bool t){
    auto mask = r.allele_mask();
    BOOST_CHECK_EQUAL(r.genotypes()[0].complex(), t);
  };
  const auto truth = vector<bool> {false, false, false, false, true, false};
  auto i = 0;
  for (const auto& rec : SingleVariantReader{"testdata/test_variants.bcf"}) {
    test_rec(rec, truth[i]);
    ++i;
  }
}

BOOST_AUTO_TEST_CASE( mixed_test ) {
  auto test_rec = [](const Variant& r, const bool t){
    auto mask = r.allele_mask();
    BOOST_CHECK_EQUAL(r.genotypes()[0].mixed(), t);
  };
  const auto truth = vector<bool> {false, false, false, false, true, false};
  auto i = 0;
  for (const auto& rec : SingleVariantReader{"testdata/test_variants.bcf"}) {
    test_rec(rec, truth[i]);
    ++i;
  }
}

BOOST_AUTO_TEST_CASE( is_variant ) {
  const auto truth = vector<bool>{true, false, true};
  const auto rec = *(SingleVariantReader{"testdata/test_variants.bcf"}.begin());
  const auto genotypes = rec.genotypes();
  for (auto i=0u; i != genotypes.size(); ++i) 
    BOOST_CHECK_EQUAL(genotypes[i].variant(), truth[i]);
}

