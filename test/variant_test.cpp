#include <boost/test/unit_test.hpp>
#include "variant_reader.h"
#include "variant.h"
#include "utils/variant_utils.h"

using namespace std;
using namespace gamgee;
using namespace boost;

BOOST_AUTO_TEST_CASE( allele_mask_simple_snp ) {
  const auto rec = (*(SingleVariantReader("testdata/test_variants.bcf").begin())); // take the first record
  const auto am = rec.allele_mask();
  BOOST_REQUIRE_EQUAL(am.size(), 2u);
  BOOST_CHECK(am[0] == AlleleType::REFERENCE);
  BOOST_CHECK(am[1] == AlleleType::SNP);
}

BOOST_AUTO_TEST_CASE( allele_mask_simple_insertion ) {
  auto it = SingleVariantReader("testdata/test_variants.bcf").begin(); 
  ++it; ++it; ++it; // take the fourth record (it is an insertion)
  const auto rec = *it;
  const auto am = rec.allele_mask();
  BOOST_REQUIRE_EQUAL(am.size(), 2u);
  BOOST_CHECK(am[0] == AlleleType::REFERENCE);
  BOOST_CHECK(am[1] == AlleleType::INSERTION);
}

BOOST_AUTO_TEST_CASE( allele_mask_simple_deletion ) {
  auto it = SingleVariantReader("testdata/test_variants.bcf").begin(); 
  ++it; ++it; // take the third record -- this would be so much simpler with a VariantBuilder....
  const auto rec = *it;
  const auto am = rec.allele_mask();
  BOOST_REQUIRE_EQUAL(am.size(), 2u);
  BOOST_CHECK(am[0] == AlleleType::REFERENCE);
  BOOST_CHECK(am[1] == AlleleType::DELETION);
}

BOOST_AUTO_TEST_CASE( allele_mask_snp_and_insertion ) {
  auto it = SingleVariantReader("testdata/test_variants.bcf").begin(); 
  ++it; ++it; ++it; ++it; // take the third record -- this would be so much simpler with a VariantBuilder....
  const auto rec = *it;
  const auto am = rec.allele_mask();
  BOOST_REQUIRE_EQUAL(am.size(), 3u);
  BOOST_CHECK(am[0] == AlleleType::REFERENCE);
  BOOST_CHECK(am[1] == AlleleType::DELETION);
  BOOST_CHECK(am[2] == AlleleType::INSERTION);
}

