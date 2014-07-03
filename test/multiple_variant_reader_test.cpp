#include "multiple_variant_reader.h"
#include "multiple_variant_iterator.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

// TODO: temp passthrough
BOOST_AUTO_TEST_CASE( multi_variant_reader )
{
  const auto truth_contigs = vector<uint32_t>{0, 1, 1, 1, 2};
  const auto truth_alignment_starts = vector<uint32_t>{10000000, 10001000, 10002000, 10003000, 10004000};
  const auto truth_n_alleles = vector<uint32_t>{2, 2, 2, 2, 3};
  MultipleVariantReader<MultipleVariantIterator> reader{{"testdata/test_variants.vcf", "testdata/test_variants.bcf"}};
  auto record_counter = 0u;
  for (const auto& vec : reader) {
    for (const auto& record : vec) {
      BOOST_CHECK_EQUAL(record.chromosome(), truth_contigs[record_counter]);
      BOOST_CHECK_EQUAL(record.alignment_start(), truth_alignment_starts[record_counter]);
      BOOST_CHECK_EQUAL(record.n_alleles(), truth_n_alleles[record_counter]);
      BOOST_CHECK_EQUAL(record.n_samples(), 3);
      BOOST_CHECK_EQUAL(record.qual(), 0);

      // testing genotype_quals accessors
      const auto gqs = record.genotype_quals();
      BOOST_CHECK_EQUAL(gqs[1][0], 35); // checking random access operators
      for_each(gqs.begin(), gqs.end(), [](const VariantFieldValue<int32_t>& x) { BOOST_CHECK_EQUAL(x[0], 35); }); // testing std library functional call at the samples level

      // testing phred_likelihoods accessors
      const auto pls = record.phred_likelihoods();
      BOOST_CHECK_EQUAL(pls[0][0], 10); // testing random access
      BOOST_CHECK_EQUAL(pls[1][1], 10);
      BOOST_CHECK_EQUAL(pls[2][2], 0);
      for(const auto& sample_pl : pls) { // testing for-each iteration at the samples level
        BOOST_CHECK_EQUAL(1, count_if(sample_pl.begin(), sample_pl.end(), [](const auto& x) { return x == 0; })); // testing std library functional call at the sample value level
        for(const auto& value_pl : sample_pl)  // testing for-each iteration at the sample value level
          BOOST_CHECK_EQUAL(value_pl, value_pl); // I just want to make sure that this level of iteration is supported, the values don't matter anymore at this point
      }

      // check genotype accessors
      BOOST_CHECK(record.is_hom_ref(1));
      BOOST_CHECK(record.is_het(0));
      BOOST_CHECK(record.is_hom_var(2));

      // test generic FLOAT
      const auto af = record.generic_float_format_field("AF");
      for_each(af.begin(), af.end(), [](const auto& s) {BOOST_CHECK_CLOSE(s[1], 2.2, 0.001);}); // checking float std library functional call at the samples level
      for (const auto& s : af) {
        for_each(s.begin(), s.end(), [](const auto& x) {BOOST_CHECK_EQUAL(x,x);}); // checking float std library functional call at the sample value level
        BOOST_CHECK_CLOSE(s[0], 3.1, 0.001);
        BOOST_CHECK_CLOSE(s[1], 2.2, 0.001);
      }

      // test generic STRING (disabled because it is not working!)
      // const auto as = record.generic_string_format_field("AS");
      // for_each(as.begin(), as.end(), [](const auto& s) {BOOST_CHECK_EQUAL(s[1], "CATE");}); // checking string std library functional call at the samples level
      // for (const auto& s : as) {
      //   for_each(s.begin(), s.end(), [](const auto& x) {BOOST_CHECK_EQUAL(x,x);}); // checking string std library functional call at the sample value level
      //   BOOST_CHECK_EQUAL(s[0], "ABA");
      //   BOOST_CHECK_EQUAL(s[1], "CATE");
      // }
      ++record_counter;
   }
  }
  std::cout << "pre counter check" << std::endl;
  BOOST_CHECK_EQUAL(record_counter, 5u);
  std::cout << "end test" << std::endl;
}
