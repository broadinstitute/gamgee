#include "variant.h"
#include "variant_reader.h"
#include "is_missing.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( single_variant_reader ) 
{
  const auto truth_contigs = vector<uint32_t>{0, 1, 1, 1, 2};
  const auto truth_alignment_starts = vector<uint32_t>{10000000, 10001000, 10002000, 10003000, 10004000};
  const auto truth_n_alleles = vector<uint32_t>{2, 2, 2, 2, 3};
  const auto truth_filter_first = vector<string>{"PASS", "PASS", "LOW_QUAL", "NOT_DEFINED", "PASS"};
  const auto truth_filter_size = vector<uint32_t>{1,1,1,1,2};
  const auto truth_quals = vector<float>{80,8.4,-1,-1,-1};
  const auto truth_ref = vector<string>{"T", "GG", "TAGTGQA", "A", "GAT"};
  const auto truth_id = vector<string>{"db2342", "rs837472", ".", ".", "."};
  const vector< vector<string>> truth_alt = {  { "C" } , {"AA"},  {"T"},  {"AGCT"},  {"G","GATAT"}};
  for (const auto& filename : {"testdata/test_variants.vcf", "testdata/test_variants.bcf"}) {
    auto record_counter = 0u;
    for (const auto& record : SingleVariantReader{filename}) {
      BOOST_CHECK_EQUAL(record.ref(), truth_ref[record_counter]);
      BOOST_CHECK_EQUAL(record.chromosome(), truth_contigs[record_counter]);
      BOOST_CHECK_EQUAL(record.alignment_start(), truth_alignment_starts[record_counter]);
      BOOST_CHECK_EQUAL(record.n_alleles(), truth_n_alleles[record_counter]);
      BOOST_CHECK_EQUAL(record.n_samples(), 3);
      BOOST_CHECK_EQUAL(record.id(), truth_id[record_counter]);
      BOOST_CHECK(record.has_filter(truth_filter_first[record_counter])); // check the has_filter member function


      // check for quals (whether missing or parsed)
      if (truth_quals[record_counter] < 0) 
        BOOST_CHECK(is_missing(record.qual()));
      else 
        BOOST_CHECK_EQUAL(record.qual(), truth_quals[record_counter]);

      auto alt = record.alt();
      BOOST_CHECK_EQUAL_COLLECTIONS(alt.begin(), alt.end(), truth_alt[record_counter].begin(), truth_alt[record_counter].end());

      // check that absent filters are really absent
      const auto absent_filter = record.has_filter("LOW_QUAL");
      if (record_counter != 2)
        BOOST_CHECK(absent_filter == false);

      // test filters getter api
      const auto filters = record.filters(); 
      BOOST_CHECK_EQUAL(filters.size(), truth_filter_size[record_counter]); // checking size member function
      BOOST_CHECK_EQUAL(filters[0], truth_filter_first[record_counter]);  // checking random access
      BOOST_CHECK_EQUAL(count_if(filters.begin(), filters.end(), [](const auto x){return x == x;}), truth_filter_size[record_counter]); // checking iteration in functional style

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
      const auto validated_expected = std::vector<bool>{record_counter == 0};
      BOOST_CHECK_EQUAL_COLLECTIONS(validated_actual.begin(), validated_actual.end(), validated_expected.begin(), validated_expected.end());

      const auto an = record.generic_integer_info_field("AN");
      BOOST_CHECK_EQUAL(an.size(), 1);
      BOOST_CHECK_EQUAL(an[0], 6);

      const auto af_actual = record.generic_float_info_field("AF");
      const auto af_expected = record_counter == 4 ? std::vector<float>{0.5, 0} : std::vector<float>{0.5};
      BOOST_CHECK_EQUAL_COLLECTIONS(af_actual.begin(), af_actual.end(), af_expected.begin(), af_expected.end());

      const auto desc_actual = record.generic_string_info_field("DESC");
      const auto desc_expected = record_counter == 0 ? std::vector<string>{"Test1,Test2"} : std::vector<string>{};
      BOOST_CHECK_EQUAL_COLLECTIONS(desc_actual.begin(), desc_actual.end(), desc_expected.begin(), desc_expected.end());

      ++record_counter;
    }
    BOOST_CHECK_EQUAL(record_counter, 5u);
  }
}

BOOST_AUTO_TEST_CASE( single_variant_reader_sites_only )  
{
  for (const auto& record : SingleVariantReader{"testdata/test_variants.vcf", vector<string>{}})  // exclude all samples (sites-only)
    BOOST_CHECK_EQUAL(record.n_samples(), 0);
}

BOOST_AUTO_TEST_CASE( missing_id_field )  
{
  const auto truth_missing = vector<bool>{false, false, true, true, true};
  auto i = 0u;
  for (const auto& record : SingleVariantReader{"testdata/test_variants.vcf"}) 
    BOOST_CHECK_EQUAL(is_missing(record.id()), truth_missing[i++]); // check that the missing and non-missing values in the vcf actually return the right missingness
  BOOST_CHECK(is_missing(""));                                      // check that an empty string is missing
  BOOST_CHECK(!is_missing("arBItrary_vaLue"));                      // check that a non-empty string is not missing
}


BOOST_AUTO_TEST_CASE( single_variant_reader_include_all_samples ) 
{
  for (const auto& record : SingleVariantReader{"testdata/test_variants.vcf", vector<string>{}, false})  // include all samples by setting include == false and passing an empty list
    BOOST_CHECK_EQUAL(record.n_samples(), 3);
}

BOOST_AUTO_TEST_CASE( single_variant_reader_including )  
{
  for (const auto& record : SingleVariantReader{"testdata/test_variants.vcf", vector<string>{"NA12878"}}) // include only NA12878
    BOOST_CHECK_EQUAL(record.n_samples(), 1);
  for (const auto& record : SingleVariantReader{"testdata/test_variants.vcf", vector<string>{"NA12878", "NA12892"}}) // include both these samples
    BOOST_CHECK_EQUAL(record.n_samples(), 2);
}

BOOST_AUTO_TEST_CASE( single_variant_reader_excluding )  
{
  for (const auto& record : SingleVariantReader{"testdata/test_variants.vcf", vector<string>{"NA12891"}, false})  // exclude only NA12891
    BOOST_CHECK_EQUAL(record.n_samples(), 2);
  for (const auto& record : SingleVariantReader{"testdata/test_variants.vcf", vector<string>{"NA12891", "NA12878"}, false})  // exclude both these samples
    BOOST_CHECK_EQUAL(record.n_samples(), 1);
}

BOOST_AUTO_TEST_CASE( single_variant_reader_missing_data )
{
  for (const auto& record : SingleVariantReader{"testdata/test_variants_missing_data.vcf", vector<string>{"0001A","0002A","0003A","0004A","0005A","0006A","0007A"}}){ // include 7 samples
    BOOST_CHECK_EQUAL(record.n_samples(), 7); 
  }
  for (const auto& record : SingleVariantReader{"testdata/test_variants_missing_data.vcf", vector<string>{"0007B", "0008A", "0009A", "0009B"}, false}){ // exclude 4 samples
    BOOST_CHECK_EQUAL(record.n_samples(), 7); 
  }
  for (const auto& record : SingleVariantReader{"testdata/test_variants_missing_data.vcf"}){ // include all 11 samples
    BOOST_CHECK_EQUAL(record.n_samples(), 11); 
  }
}
