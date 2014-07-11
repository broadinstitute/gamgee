#include "multiple_variant_reader.h"
#include "multiple_variant_iterator.h"
#include "is_missing.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( multi_variant_reader_validation )
{
  const std::vector<std::string> filenames1{"testdata/test_variants.vcf", "testdata/extra_header.vcf"};
  const std::vector<std::string> filenames2{"testdata/test_variants.vcf", "testdata/missing_header.vcf"};

  // validate mismatched headers by default
  for (const auto filenames_v : {filenames1, filenames2})
    BOOST_CHECK_THROW(
      auto reader = MultipleVariantReader<MultipleVariantIterator>(filenames_v),
      std::runtime_error
      // this does not work for some reason
      // MultipleVariantReader::HeaderException
    );

  // don't validate mismatched headers
  for (const auto filenames_v : {filenames1, filenames2})
    auto reader = MultipleVariantReader<MultipleVariantIterator>(filenames_v, false);
}

// copied from variant_reader_test
BOOST_AUTO_TEST_CASE( multi_variant_reader )
{
  const auto truth_contigs = vector<uint32_t>{0, 1, 1, 1, 2};
  const auto truth_alignment_starts = vector<uint32_t>{10000000, 10001000, 10002000, 10003000, 10004000};
  const auto truth_n_alleles = vector<uint32_t>{2, 2, 2, 2, 3};
  const auto truth_filter_first = vector<string>{"PASS", "PASS", "LOW_QUAL", "NOT_DEFINED", "PASS"};
  const auto truth_filter_size = vector<uint32_t>{1,1,1,1,2};
  const auto truth_quals = vector<float>{80,8.4,-1,-1,-1};
  const auto truth_ref = vector<string>{"T", "GG", "TAGTGQA", "A", "GAT"};
  const vector< vector<string>> truth_alt = {  { "C" } , {"AA"},  {"T"},  {"AGCT"},  {"G","GATAT"}};
  MultipleVariantReader<MultipleVariantIterator> reader{{"testdata/test_variants.vcf", "testdata/test_variants.bcf"}, false};
  auto position_counter = 0u;
  for (const auto& vec : reader) {
    for (const auto& record : vec) {
      BOOST_CHECK_EQUAL(record.ref(), truth_ref[position_counter]);
      BOOST_CHECK_EQUAL(record.chromosome(), truth_contigs[position_counter]);
      BOOST_CHECK_EQUAL(record.alignment_start(), truth_alignment_starts[position_counter]);
      BOOST_CHECK_EQUAL(record.n_alleles(), truth_n_alleles[position_counter]);
      BOOST_CHECK_EQUAL(record.n_samples(), 3);
      BOOST_CHECK(record.has_filter(truth_filter_first[position_counter])); // check the has_filter member function


      // check for quals (whether missing or parsed)
      if (truth_quals[position_counter] < 0)
	BOOST_CHECK(is_missing(record.qual()));
      else
	BOOST_CHECK_EQUAL(record.qual(), truth_quals[position_counter]);

      auto alt = record.alt();
      BOOST_CHECK_EQUAL_COLLECTIONS(alt.begin(), alt.end(), truth_alt[position_counter].begin(), truth_alt[position_counter].end());

      // check that absent filters are really absent
      const auto absent_filter = record.has_filter("LOW_QUAL");
      if (position_counter != 2)
	BOOST_CHECK(absent_filter == false);

      // test filters getter api
      const auto filters = record.filters();
      BOOST_CHECK_EQUAL(filters.size(), truth_filter_size[position_counter]); // checking size member function
      BOOST_CHECK_EQUAL(filters[0], truth_filter_first[position_counter]);  // checking random access
      BOOST_CHECK_EQUAL(count_if(filters.begin(), filters.end(), [](const auto x){return x == x;}), truth_filter_size[position_counter]); // checking iteration in functional style

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

      const auto validated_actual = record.generic_boolean_info_field("VALIDATED");
      const auto validated_expected = std::vector<bool>{position_counter == 0};
      BOOST_CHECK_EQUAL_COLLECTIONS(validated_actual.begin(), validated_actual.end(), validated_expected.begin(), validated_expected.end());

      const auto an = record.generic_integer_info_field("AN");
      BOOST_CHECK_EQUAL(an.size(), 1);
      BOOST_CHECK_EQUAL(an[0], 6);

      const auto af_actual = record.generic_float_info_field("AF");
      const auto af_expected = position_counter == 4 ? std::vector<float>{0.5, 0} : std::vector<float>{0.5};
      BOOST_CHECK_EQUAL_COLLECTIONS(af_actual.begin(), af_actual.end(), af_expected.begin(), af_expected.end());

      const auto desc_actual = record.generic_string_info_field("DESC");
      const auto desc_expected = position_counter == 0 ? std::vector<string>{"Test1,Test2"} : std::vector<string>{};
      BOOST_CHECK_EQUAL_COLLECTIONS(desc_actual.begin(), desc_actual.end(), desc_expected.begin(), desc_expected.end());
    }
    ++position_counter;
  }
  BOOST_CHECK_EQUAL(position_counter, 5u);
}
