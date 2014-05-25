#include "variant_reader.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( single_variant_reader ) 
{
  const auto truth_contigs = vector<uint32_t>{0, 1, 1, 1, 2};
  const auto truth_alignment_starts = vector<uint32_t>{10000000, 10001000, 10002000, 10003000, 10004000};
  const auto truth_n_alleles = vector<uint32_t>{2, 2, 2, 2, 3};
  for (const auto& filename : {"testdata/test_variants.vcf"}) {
    auto record_counter = 0u;
    for (const auto& record : SingleVariantReader{filename}) {
      BOOST_CHECK_EQUAL(record.chromosome(), truth_contigs[record_counter]);
      BOOST_CHECK_EQUAL(record.alignment_start(), truth_alignment_starts[record_counter]);
      BOOST_CHECK_EQUAL(record.n_alleles(), truth_n_alleles[record_counter]);
      BOOST_CHECK_EQUAL(record.n_samples(), 3);
      BOOST_CHECK_EQUAL(record.qual(), 0);
      const auto gqs = record.genotype_quals();
      BOOST_CHECK_EQUAL(gqs[1][0], 35); // checking random access operators
      for_each(gqs.begin(), gqs.end(), [](const VariantFieldValue<int32_t>& x) { BOOST_CHECK_EQUAL(x[0], 35); });
      const auto pls = record.phred_likelihoods();
      for(const auto sample_pl : pls) 
        BOOST_CHECK_EQUAL(sample_pl[1], 0);
      ++record_counter;
    }
    BOOST_CHECK_EQUAL(record_counter, 5u);
  }
}

