#include "variant.h"
#include "variant_reader.h"
#include "multiple_variant_reader.h"
#include "multiple_variant_iterator.h"
#include "indexed_variant_reader.h"
#include "indexed_variant_iterator.h"
#include "missing.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

constexpr auto FLOAT_COMPARISON_THRESHOLD = 0.0001f;

const auto truth_chromosome        = vector<uint32_t>{0, 1, 1, 1, 2};
const auto truth_alignment_starts  = vector<uint32_t>{10000000, 10001000, 10002000, 10003000, 10004000};
const auto truth_alignment_stops   = vector<uint32_t>{10000000, 10001001, 10002006, 10003000, 10004002};
const auto truth_n_alleles         = vector<uint32_t>{2, 2, 2, 2, 3};
const auto truth_filter_name       = vector<string>{"PASS", "PASS", "LOW_QUAL", "NOT_DEFINED", "PASS"};
const auto truth_filter_size       = vector<uint32_t>{1,1,1,1,2};
const auto truth_quals             = vector<float>{80,8.4,-1,-1,-1};
const auto truth_ref               = vector<string>{"T", "GG", "TAGTGQA", "A", "GAT"};
const auto truth_alt               = vector< vector<string>> {  { "C" } , {"AA"},  {"T"},  {"AGCT"},  {"G","GATAT"}};
const auto truth_high_quality_hets = boost::dynamic_bitset<>{std::string{"001"}};
const auto truth_id                = vector<string>{"db2342", "rs837472", ".", ".", "."};
const auto truth_shared_af         = vector<vector<float>>{{0.5}, {0.5}, {0.5}, {0.5}, {0.5, 0.0}};
const auto truth_shared_an         = vector<vector<int32_t>>{{6}, {6}, {6}, {6}, {6}};
const auto truth_shared_desc       = vector<vector<string>>{{"Test1,Test2"}, {}, {}, {}, {}};
const auto truth_shared_validated  = vector<bool>{true, false, true, false, false};
const auto truth_gq                = vector<vector<uint32_t>>{{25,12,650}, {35,35,35}, {35,35,35}, {35,35,35}, {35,35,35}};
const auto truth_af                = vector<float> { 3.1,2.2 };
const auto truth_pl                = vector<vector<vector<uint32_t>>>{
  {{10,0,100      }, {0,10,1000      }, {10,100,0}      },
  {{10,0,100      }, {0,10,100       }, {10,100,0}      },
  {{10,0,100      }, {0,10,2000000000}, {10,100,0}      },
  {{10,0,100      }, {0,10,100       }, {10,100,0}      },
  {{10,0,100,2,4,8}, {0,10,100,2,4,8 }, {10,100,0,2,4,8}}
};
const auto truth_as                = vector<vector<string>>{ 
  {"ABA","CA","XISPAFORRO"}, 
  {"ABA","ABA","ABA"}, 
  {"ABA","ABA","ABA"}, 
  {"ABA","ABA","ABA"}, 
  {"ABA","ABA","ABA"} 
};




boost::dynamic_bitset<> high_qual_hets(const Variant& record) {  // filter all hets that have GQ > 20
  const auto genotypes = record.genotypes(); // a "vector-like" with the genotypes of all samples in this record
  const auto gqs = record.integer_individual_field("GQ"); // a "vector-like" with all the GQs of all samples in this record
  const auto hets = Variant::select_if(genotypes.begin(), genotypes.end(), [](const auto& g) { return g.het(); }); // returns a bit set with all hets marked with 1's
  const auto pass_gqs = Variant::select_if(gqs.begin(), gqs.end(), [](const auto& gq) { return gq[0] > 20; }); // returns a bit set with every sample with gq > 20 marked with 1's
  return hets & pass_gqs; // returns a bit set with all the samples that are het and have gq > 20
}

void check_variant_basic_api(const Variant& record, const uint32_t truth_index) {
  BOOST_CHECK_EQUAL(record.ref(), truth_ref[truth_index]);
  BOOST_CHECK_EQUAL(record.chromosome(), truth_chromosome[truth_index]);
  BOOST_CHECK_EQUAL(record.alignment_start(), truth_alignment_starts[truth_index]);
  BOOST_CHECK_EQUAL(record.alignment_stop(), truth_alignment_stops[truth_index]);
  BOOST_CHECK_EQUAL(record.n_alleles(), truth_n_alleles[truth_index]);
  BOOST_CHECK_EQUAL(record.n_samples(), 3);
  BOOST_CHECK_EQUAL(record.id(), truth_id[truth_index]);
}

void check_quals_api(const Variant& record, const uint32_t truth_index) {
  if (truth_quals[truth_index] < 0) 
    BOOST_CHECK(missing(record.qual()));
  else 
    BOOST_CHECK_EQUAL(record.qual(), truth_quals[truth_index]); 
}

void check_alt_api(const Variant& record, const uint32_t truth_index) {
  auto alt = record.alt();
  BOOST_CHECK_EQUAL_COLLECTIONS(alt.begin(), alt.end(), truth_alt[truth_index].begin(), truth_alt[truth_index].end());
}

void check_filters_api(const Variant& record, const uint32_t truth_index) {
  BOOST_CHECK(record.has_filter(truth_filter_name[truth_index])); // check the has_filter member function
  const auto low_qual_filter_present = record.has_filter("LOW_QUAL");
  if (truth_index == 2)
    BOOST_CHECK(low_qual_filter_present == true);
  else
    BOOST_CHECK(low_qual_filter_present == false);
  const auto filters = record.filters(); 
  BOOST_CHECK_EQUAL(filters.size(), truth_filter_size[truth_index]); // checking size member function
  BOOST_CHECK_EQUAL(filters[0], truth_filter_name[truth_index]);  // checking random access
  BOOST_CHECK_EQUAL(count_if(filters.begin(), filters.end(), [](const auto x){return x == x;}), truth_filter_size[truth_index]); // checking iteration in functional style
}

void check_genotype_quals_api(const Variant& record, const uint32_t truth_index) {
  const auto gqs = record.integer_individual_field("GQ");
  // direct access
  for (auto i = 0u; i != record.n_samples(); ++i) 
    BOOST_CHECK_EQUAL(gqs[i][0], truth_gq[truth_index][i]);
  // functional style access
  auto i = 0;
  for_each(gqs.begin(), gqs.end(), [&truth_index, &i](const IndividualFieldValue<int32_t>& x) { BOOST_CHECK_EQUAL(x[0], truth_gq[truth_index][i++]); }); // testing std library functional call at the samples level
}

void check_phred_likelihoods_api(const Variant& record, const uint32_t truth_index) {
  const auto pls = record.integer_individual_field("PL");
  // direct access
  for(auto i = 0u; i != record.n_samples(); ++i) 
    for (auto j = 0u; j != pls[i].size(); ++j) 
      BOOST_CHECK_EQUAL(pls[i][j], truth_pl[truth_index][i][j]);
  // functional style access
  for(const auto& sample_pl : pls) { // testing for-each iteration at the samples level
    BOOST_CHECK_EQUAL(1, count_if(sample_pl.begin(), sample_pl.end(), [](const auto& x) { return x == 0; })); // testing std library functional call at the sample value level
    for(const auto& value_pl : sample_pl)  // testing for-each iteration at the sample value level
      BOOST_CHECK_EQUAL(value_pl, value_pl); // I just want to make sure that this level of iteration is supported, the values don't matter anymore at this point
  }
}      

void check_individual_field_api(const Variant& record, const uint32_t truth_index) {
  const auto gq_int    = record.individual_field_as_integer("GQ");
  const auto gq_float  = record.individual_field_as_float("GQ");
  const auto gq_string = record.individual_field_as_string("GQ");
  const auto af_int    = record.individual_field_as_integer("AF");
  const auto af_float  = record.individual_field_as_float("AF");
  const auto af_string = record.individual_field_as_string("AF");
  const auto pl_int    = record.individual_field_as_integer("PL");
  const auto pl_float  = record.individual_field_as_float("PL");
  const auto pl_string = record.individual_field_as_string("PL"); 
  const auto as_int    = record.individual_field_as_integer("AS"); // this is a string field, we should be able to create the object but not access it's elements due to lazy initialization
  const auto as_float  = record.individual_field_as_float("AS");   // this is a string field, we should be able to create the object but not access it's elements due to lazy initialization
  const auto as_string = record.individual_field_as_string("AS");
  
  // index based API
  const auto header = record.header();
  const auto gq_idx = header.field_index("GQ");
  const auto af_idx = header.field_index("AF");
  const auto pl_idx = header.field_index("PL");
  const auto as_idx = header.field_index("AS");
  const auto gq_int_idx    = record.individual_field_as_integer(gq_idx);
  const auto gq_float_idx  = record.individual_field_as_float(gq_idx);
  const auto gq_string_idx = record.individual_field_as_string(gq_idx);
  const auto af_int_idx    = record.individual_field_as_integer(af_idx);
  const auto af_float_idx  = record.individual_field_as_float(af_idx);
  const auto af_string_idx = record.individual_field_as_string(af_idx);
  const auto pl_int_idx    = record.individual_field_as_integer(pl_idx);
  const auto pl_float_idx  = record.individual_field_as_float(pl_idx);
  const auto pl_string_idx = record.individual_field_as_string(pl_idx); 
  const auto as_int_idx    = record.individual_field_as_integer(as_idx); // this is a string field, we should be able to create the object but not access it's elements due to lazy initialization
  const auto as_float_idx  = record.individual_field_as_float(as_idx);   // this is a string field, we should be able to create the object but not access it's elements due to lazy initialization
  const auto as_string_idx = record.individual_field_as_string(as_idx);


  // test the conversions using unforgiving API
  BOOST_CHECK_THROW(record.float_individual_field("GQ"), runtime_error);
  BOOST_CHECK_THROW(record.string_individual_field("GQ"), runtime_error);
  BOOST_CHECK_THROW(record.integer_individual_field("AF"), runtime_error);
  BOOST_CHECK_THROW(record.string_individual_field("AF"), runtime_error);
  BOOST_CHECK_THROW(record.float_individual_field("PL"), runtime_error);
  BOOST_CHECK_THROW(record.string_individual_field("PL"), runtime_error);
  BOOST_CHECK_THROW(record.integer_individual_field("AS"), runtime_error);
  BOOST_CHECK_THROW(record.float_individual_field("AS"), runtime_error);
  BOOST_CHECK_THROW(record.float_individual_field(gq_idx), runtime_error);
  BOOST_CHECK_THROW(record.string_individual_field(gq_idx), runtime_error);
  BOOST_CHECK_THROW(record.integer_individual_field(af_idx), runtime_error);
  BOOST_CHECK_THROW(record.string_individual_field(af_idx), runtime_error);
  BOOST_CHECK_THROW(record.float_individual_field(pl_idx), runtime_error);
  BOOST_CHECK_THROW(record.string_individual_field(pl_idx), runtime_error);
  BOOST_CHECK_THROW(record.integer_individual_field(as_idx), runtime_error);
  BOOST_CHECK_THROW(record.float_individual_field(as_idx), runtime_error);
  
  BOOST_CHECK(gq_int == record.integer_individual_field("GQ"));
  BOOST_CHECK(af_float == record.float_individual_field("AF"));
  BOOST_CHECK(pl_int == record.integer_individual_field("PL"));
  BOOST_CHECK(as_string == record.string_individual_field("AS"));
  BOOST_CHECK(gq_int != record.individual_field_as_integer("AF"));
  BOOST_CHECK(pl_float != record.individual_field_as_float("GQ"));
  BOOST_CHECK(as_string != record.individual_field_as_string("PL"));

  for(auto i=0u; i != record.n_samples(); ++i) {
    BOOST_CHECK_EQUAL(gq_int[i][0], truth_gq[truth_index][i]);
    BOOST_CHECK_EQUAL(gq_int_idx[i][0], truth_gq[truth_index][i]);
    BOOST_CHECK_CLOSE(gq_float[i][0], float(truth_gq[truth_index][i]), FLOAT_COMPARISON_THRESHOLD);
    BOOST_CHECK_CLOSE(gq_float_idx[i][0], float(truth_gq[truth_index][i]), FLOAT_COMPARISON_THRESHOLD);
    BOOST_CHECK_EQUAL(gq_string[i][0], to_string(truth_gq[truth_index][i]));
    BOOST_CHECK_EQUAL(gq_string_idx[i][0], to_string(truth_gq[truth_index][i]));
    BOOST_REQUIRE_EQUAL(af_int[i].size(), truth_af.size()); // require otherwise next line may segfault
    BOOST_REQUIRE_EQUAL(af_int_idx[i].size(), truth_af.size()); // require otherwise next line may segfault
    for (auto j=0u; j!= af_int[i].size(); ++j) {
      BOOST_CHECK_EQUAL(af_int[i][j], int32_t(truth_af[j]));
      BOOST_CHECK_EQUAL(af_int_idx[i][j], int32_t(truth_af[j]));
      BOOST_CHECK_CLOSE(af_float[i][j], truth_af[j], FLOAT_COMPARISON_THRESHOLD);
      BOOST_CHECK_CLOSE(af_float_idx[i][j], truth_af[j], FLOAT_COMPARISON_THRESHOLD);
      BOOST_CHECK_EQUAL(af_string[i][j], to_string(truth_af[j]));
      BOOST_CHECK_EQUAL(af_string_idx[i][j], to_string(truth_af[j]));
    }
    BOOST_REQUIRE_EQUAL(pl_int[i].size(), truth_pl[truth_index][i].size()); // require otherwise next line may segfault
    BOOST_REQUIRE_EQUAL(pl_int_idx[i].size(), truth_pl[truth_index][i].size()); // require otherwise next line may segfault
    for (auto j=0u; j!= pl_int[i].size(); ++j) { 
      BOOST_CHECK_EQUAL(pl_int[i][j], truth_pl[truth_index][i][j]);
      BOOST_CHECK_EQUAL(pl_int_idx[i][j], truth_pl[truth_index][i][j]);
      BOOST_CHECK_CLOSE(pl_float[i][j], float(truth_pl[truth_index][i][j]), FLOAT_COMPARISON_THRESHOLD);
      BOOST_CHECK_CLOSE(pl_float_idx[i][j], float(truth_pl[truth_index][i][j]), FLOAT_COMPARISON_THRESHOLD);
      BOOST_CHECK_EQUAL(pl_string[i][j], to_string(truth_pl[truth_index][i][j]));
      BOOST_CHECK_EQUAL(pl_string_idx[i][j], to_string(truth_pl[truth_index][i][j]));
    }
    BOOST_CHECK_EQUAL(as_string[i][0], truth_as[truth_index][i]);
    BOOST_CHECK_EQUAL(as_string_idx[i][0], truth_as[truth_index][i]);
  }
  BOOST_CHECK_THROW(as_float[0][0], invalid_argument); 
  BOOST_CHECK_THROW(as_float_idx[0][0], invalid_argument); 
  BOOST_CHECK_THROW(as_int[0][0], invalid_argument); 
  BOOST_CHECK_THROW(as_int_idx[0][0], invalid_argument); 
  BOOST_CHECK(missing(record.integer_individual_field("NON_EXISTING")));
  BOOST_CHECK(missing(record.integer_individual_field(-1)));
  BOOST_CHECK(missing(record.float_individual_field("NON_EXISTING")));
  BOOST_CHECK(missing(record.float_individual_field(-1)));
  BOOST_CHECK(missing(record.string_individual_field("NON_EXISTING")));
  BOOST_CHECK(missing(record.string_individual_field(-1)));
  BOOST_CHECK_THROW(record.float_individual_field("NON_EXISTING")[0], out_of_range);
  BOOST_CHECK_THROW(record.float_individual_field(-1)[0], out_of_range);
}

void check_shared_field_api(const Variant& record, const uint32_t truth_index) {
  const auto header = record.header();
  BOOST_CHECK_EQUAL(record.boolean_shared_field("VALIDATED"), truth_shared_validated[truth_index]);
  const auto an = record.integer_shared_field("AN"); // test the string based api
  BOOST_CHECK_EQUAL_COLLECTIONS(an.begin(), an.end(), truth_shared_an[truth_index].begin(), truth_shared_an[truth_index].end());
  const auto an_idx = record.integer_shared_field(header.field_index("AN")); // test the index based api
  BOOST_CHECK_EQUAL_COLLECTIONS(an_idx.begin(), an_idx.end(), truth_shared_an[truth_index].begin(), truth_shared_an[truth_index].end());
  const auto af = record.float_shared_field("AF");
  BOOST_CHECK_EQUAL_COLLECTIONS(af.begin(), af.end(), truth_shared_af[truth_index].begin(), truth_shared_af[truth_index].end());
  const auto af_idx = record.float_shared_field(header.field_index("AF")); // test the index based api
  BOOST_CHECK_EQUAL_COLLECTIONS(af_idx.begin(), af_idx.end(), truth_shared_af[truth_index].begin(), truth_shared_af[truth_index].end());
  const auto desc = record.string_shared_field("DESC");
  BOOST_CHECK_EQUAL_COLLECTIONS(desc.begin(), desc.end(), truth_shared_desc[truth_index].begin(), truth_shared_desc[truth_index].end());
  const auto desc_idx = record.string_shared_field(header.field_index("DESC")); // test the index based api
  BOOST_CHECK_EQUAL_COLLECTIONS(desc_idx.begin(), desc_idx.end(), truth_shared_desc[truth_index].begin(), truth_shared_desc[truth_index].end());
  // check non-existing missing values
  BOOST_CHECK(missing(record.boolean_shared_field("NON_EXISTING")));
  BOOST_CHECK(missing(record.integer_shared_field("NON_EXISTING")));
  BOOST_CHECK(missing(record.float_shared_field("NON_EXISTING")));  
  BOOST_CHECK(missing(record.string_shared_field("NON_EXISTING"))); 
  BOOST_CHECK(missing(record.boolean_shared_field(-1)));
  BOOST_CHECK(missing(record.integer_shared_field(-1)));
  BOOST_CHECK(missing(record.float_shared_field(-1)));  
  BOOST_CHECK(missing(record.string_shared_field(-1))); 
  // check type conversions in the unforgiving API
  BOOST_CHECK_THROW(record.float_shared_field("AN"), runtime_error);
  BOOST_CHECK_THROW(record.string_shared_field("AN"), runtime_error);
  BOOST_CHECK_THROW(record.integer_shared_field("AF"), runtime_error);
  BOOST_CHECK_THROW(record.string_shared_field("AF"), runtime_error);
  BOOST_CHECK_THROW(record.integer_shared_field("DESC"), runtime_error);
  BOOST_CHECK_THROW(record.float_shared_field("DESC"), runtime_error);
  BOOST_CHECK_THROW(record.float_shared_field(header.field_index("AN")), runtime_error);
  BOOST_CHECK_THROW(record.string_shared_field(header.field_index("AN")), runtime_error);
  BOOST_CHECK_THROW(record.integer_shared_field(header.field_index("AF")), runtime_error);
  BOOST_CHECK_THROW(record.string_shared_field(header.field_index("AF")), runtime_error);
  BOOST_CHECK_THROW(record.integer_shared_field(header.field_index("DESC")), runtime_error);
  BOOST_CHECK_THROW(record.float_shared_field(header.field_index("DESC")), runtime_error);
  // check conversions on the nice API
  const auto an_float = record.shared_field_as_float("AN");
  const auto an_string = record.shared_field_as_string("AN");
  for (auto i=0u; i != an_float.size(); ++i) {
    BOOST_CHECK_EQUAL(an_float[i], float(truth_shared_an[truth_index][i]));
    BOOST_CHECK_EQUAL(an_string[i], to_string(truth_shared_an[truth_index][i]));
  }
  const auto an_float_idx = record.shared_field_as_float(header.field_index("AN"));
  const auto an_string_idx = record.shared_field_as_string(header.field_index("AN"));
  for (auto i=0u; i != an_float.size(); ++i) {
    BOOST_CHECK_EQUAL(an_float_idx[i], float(truth_shared_an[truth_index][i]));
    BOOST_CHECK_EQUAL(an_string_idx[i], to_string(truth_shared_an[truth_index][i]));
  }
  const auto af_integer = record.shared_field_as_integer("AF");
  const auto af_string = record.shared_field_as_string("AF");
  for (auto i=0u; i != af_integer.size(); ++i) {
    BOOST_CHECK_EQUAL(af_integer[i], int32_t(truth_shared_af[truth_index][i]));
    BOOST_CHECK_EQUAL(af_string[i], to_string(truth_shared_af[truth_index][i]));
  }
  const auto af_integer_idx = record.shared_field_as_integer(header.field_index("AF"));
  const auto af_string_idx = record.shared_field_as_string(header.field_index("AF"));
  for (auto i=0u; i != af_integer.size(); ++i) {
    BOOST_CHECK_EQUAL(af_integer_idx[i], int32_t(truth_shared_af[truth_index][i]));
    BOOST_CHECK_EQUAL(af_string_idx[i], to_string(truth_shared_af[truth_index][i]));
  }
  const auto desc_index = header.field_index("DESC");
  const auto desc_bool = record.boolean_shared_field("DESC");
  const auto desc_integer = record.shared_field_as_integer("DESC");
  const auto desc_float = record.shared_field_as_float("DESC");
  const auto desc_bool_idx = record.boolean_shared_field(desc_index);
  const auto desc_integer_idx = record.shared_field_as_integer(desc_index);
  const auto desc_float_idx = record.shared_field_as_float(desc_index);
  if (truth_shared_desc[truth_index].empty()) {
    BOOST_CHECK(!desc_bool);
    BOOST_CHECK(missing(desc_integer));
    BOOST_CHECK(missing(desc_float));
    BOOST_CHECK(!desc_bool_idx);
    BOOST_CHECK(missing(desc_integer_idx));
    BOOST_CHECK(missing(desc_float_idx));
    BOOST_CHECK(missing(record.shared_field_as_string("DESC"))); // check that an existing tag in the header can be missing 
    BOOST_CHECK(missing(record.shared_field_as_string(desc_index))); // check that an existing tag in the header can be missing 
    BOOST_CHECK_THROW(desc_float[0], out_of_range);
  }
  else {
    BOOST_CHECK(desc_bool);
    BOOST_CHECK_THROW(desc_float[0], invalid_argument);
    BOOST_CHECK_THROW(desc_integer[0], invalid_argument);
    BOOST_CHECK(desc_bool_idx);
    BOOST_CHECK_THROW(desc_float_idx[0], invalid_argument);
    BOOST_CHECK_THROW(desc_integer_idx[0], invalid_argument);
    BOOST_CHECK(desc == record.shared_field_as_string("DESC"));  
    BOOST_CHECK(desc_idx == record.shared_field_as_string(desc_index));  
  }
  BOOST_CHECK(an == record.shared_field_as_integer("AN")); 
  BOOST_CHECK(an == record.shared_field_as_integer(header.field_index("AN"))); 
  BOOST_CHECK(af == record.shared_field_as_float("AF")); 
  BOOST_CHECK(af == record.shared_field_as_float(header.field_index("AF"))); 
  BOOST_CHECK(an != record.shared_field_as_integer("AF")); 
  BOOST_CHECK(af != record.shared_field_as_float(header.field_index("AN"))); 
}

void check_genotype_api(const Variant& record, const uint32_t truth_index) {
  BOOST_CHECK_EQUAL(high_qual_hets(record), truth_high_quality_hets);
  const auto gt_for_all_samples = record.genotypes();
  BOOST_CHECK(gt_for_all_samples[1].hom_ref());
  BOOST_CHECK(gt_for_all_samples[0].het());
  BOOST_CHECK(gt_for_all_samples[2].hom_var());
  for (const auto& gt_for_single_sample: gt_for_all_samples) {
    BOOST_CHECK_EQUAL(gt_for_single_sample.size(), 2);
    const auto alleles = gt_for_single_sample.alleles_strings();
    const auto allele_keys = gt_for_single_sample.alleles_keys();
    BOOST_CHECK_EQUAL(gt_for_single_sample.allele_string(0), alleles[0]);
    BOOST_CHECK_EQUAL(gt_for_single_sample.allele_string(1), alleles[1]);
    BOOST_CHECK_EQUAL(gt_for_single_sample.allele_key(0), allele_keys[0]);
    BOOST_CHECK_EQUAL(gt_for_single_sample.allele_key(1), allele_keys[1]);
    BOOST_CHECK_EQUAL(gt_for_single_sample[0], allele_keys[0]);
    BOOST_CHECK_EQUAL(gt_for_single_sample[1], allele_keys[1]);
    BOOST_CHECK(!missing(gt_for_single_sample));
    BOOST_CHECK_EQUAL(gt_for_single_sample.size(), 2);
    if (gt_for_single_sample.het()) {
      BOOST_CHECK(!gt_for_single_sample.hom_ref());
      BOOST_CHECK(!gt_for_single_sample.hom_var());
      BOOST_CHECK_NE(allele_keys[0], allele_keys[1]);
      BOOST_CHECK_EQUAL(allele_keys[0], 0);
      BOOST_CHECK_EQUAL(allele_keys[1], 1);
      BOOST_CHECK_NE(alleles[0], alleles[1]);
      BOOST_CHECK_EQUAL(alleles[0], truth_ref[truth_index]);
      BOOST_CHECK_EQUAL(alleles[1], truth_alt[truth_index][0]);
      BOOST_CHECK_EQUAL(gt_for_single_sample.fast_diploid_key_generation(), 0x00000001);
    }
    if (gt_for_single_sample.hom_ref()) {
      BOOST_CHECK(!gt_for_single_sample.het());
      BOOST_CHECK(!gt_for_single_sample.hom_var());
      BOOST_CHECK_EQUAL(allele_keys[0], allele_keys[1]);
      BOOST_CHECK_EQUAL(allele_keys[0], 0);
      BOOST_CHECK_EQUAL(allele_keys[1], 0);
      BOOST_CHECK_EQUAL(alleles[0], alleles[1]);
      BOOST_CHECK_EQUAL(alleles[0], truth_ref[truth_index]);
      BOOST_CHECK_EQUAL(alleles[1], truth_ref[truth_index]);
      BOOST_CHECK_EQUAL(gt_for_single_sample.fast_diploid_key_generation(), 0x00000000);
    }
    if (gt_for_single_sample.hom_var()) {
      BOOST_CHECK(!gt_for_single_sample.het());
      BOOST_CHECK(!gt_for_single_sample.hom_ref());
      BOOST_CHECK_EQUAL(allele_keys[0], allele_keys[1]);
      BOOST_CHECK_EQUAL(allele_keys[0], 1);
      BOOST_CHECK_EQUAL(allele_keys[1], 1);
      BOOST_CHECK_EQUAL(alleles[0], alleles[1]);
      BOOST_CHECK_EQUAL(alleles[0], truth_alt[truth_index][0]);
      BOOST_CHECK_EQUAL(alleles[1], truth_alt[truth_index][0]);
      BOOST_CHECK_EQUAL(gt_for_single_sample.fast_diploid_key_generation(), 0x00010001);
    }
    if (truth_index == 0) {
      BOOST_CHECK(gt_for_all_samples[0] == gt_for_all_samples[0]);
      BOOST_CHECK(!(gt_for_all_samples[0] != gt_for_all_samples[0]));
      BOOST_CHECK(gt_for_all_samples[0] != gt_for_all_samples[1]);
      BOOST_CHECK(!(gt_for_all_samples[0] == gt_for_all_samples[1]));
      BOOST_CHECK(gt_for_all_samples[0][1] == gt_for_all_samples[2][0]);
      BOOST_CHECK(!(gt_for_all_samples[0][1] != gt_for_all_samples[2][0]));
      BOOST_CHECK(gt_for_all_samples[0][1] == gt_for_all_samples[0][1]);
      BOOST_CHECK(!(gt_for_all_samples[0][1] != gt_for_all_samples[0][1]));
      BOOST_CHECK(gt_for_all_samples[0][1] != gt_for_all_samples[1][0]);
      BOOST_CHECK(!(gt_for_all_samples[0][1] == gt_for_all_samples[1][0]));
    }
  }
}

void generic_variant_reader_test(const std::function<void(const Variant&, const uint32_t)>& fun) {
  for (const auto& filename : {"testdata/test_variants.vcf", "testdata/test_variants.bcf"}) {
    auto truth_index = 0u;
    for (const auto& record : SingleVariantReader{filename}) {
      fun(record, truth_index);
      ++truth_index;
    }
  }
}

BOOST_AUTO_TEST_CASE( basic_api )             { generic_variant_reader_test(check_variant_basic_api);     }
BOOST_AUTO_TEST_CASE( quals_api )             { generic_variant_reader_test(check_quals_api);             }
BOOST_AUTO_TEST_CASE( alt_api )               { generic_variant_reader_test(check_alt_api);               }
BOOST_AUTO_TEST_CASE( filters_api )           { generic_variant_reader_test(check_filters_api);           }
BOOST_AUTO_TEST_CASE( genotype_quals_api )    { generic_variant_reader_test(check_genotype_quals_api);    }
BOOST_AUTO_TEST_CASE( phred_likelihoods_api ) { generic_variant_reader_test(check_phred_likelihoods_api); }
BOOST_AUTO_TEST_CASE( individual_fields_api ) { generic_variant_reader_test(check_individual_field_api);  }
BOOST_AUTO_TEST_CASE( shared_fields_api )     { generic_variant_reader_test(check_shared_field_api);      }
BOOST_AUTO_TEST_CASE( genotype_api )          { generic_variant_reader_test(check_genotype_api);          }

BOOST_AUTO_TEST_CASE( missing_id_field )
{
  const auto truth_missing = vector<bool>{false, false, true, true, true};
  auto i = 0u;
  for (const auto& record : SingleVariantReader{"testdata/test_variants.vcf"})
    BOOST_CHECK_EQUAL(missing(record.id()), truth_missing[i++]); // check that the missing and non-missing values in the vcf actually return the right missingness
  BOOST_CHECK(missing(""));                                      // check that an empty string is missing
  BOOST_CHECK(!missing("arBItrary_vaLue"));                      // check that a non-empty string is not missing
}

void single_variant_reader_sample_test(const string filename, const vector<string> samples, const bool include, const int desired_samples) {
  for (const auto& record : SingleVariantReader{filename, samples, include})
    BOOST_CHECK_EQUAL(record.n_samples(), desired_samples);
}

BOOST_AUTO_TEST_CASE( single_variant_reader_sites_only )  
{
  single_variant_reader_sample_test("testdata/test_variants.vcf", vector<string>{}, true, 0); // exclude all samples (sites-only)
}

BOOST_AUTO_TEST_CASE( single_variant_reader_include_all_samples ) 
{
  single_variant_reader_sample_test("testdata/test_variants.vcf", vector<string>{}, false, 3);  // include all samples by setting include == false and passing an empty list
}

BOOST_AUTO_TEST_CASE( single_variant_reader_including )  
{
  single_variant_reader_sample_test("testdata/test_variants.vcf", vector<string>{"NA12878"}, true, 1); // include only NA12878
  single_variant_reader_sample_test("testdata/test_variants.vcf", vector<string>{"NA12878", "NA12892"}, true, 2); // include both these samples
}

BOOST_AUTO_TEST_CASE( single_variant_reader_excluding )  
{
  single_variant_reader_sample_test("testdata/test_variants.vcf", vector<string>{"NA12891"}, false, 2);  // exclude only NA12891
  single_variant_reader_sample_test("testdata/test_variants.vcf", vector<string>{"NA12891", "NA12878"}, false, 1);  // exclude both these samples
}

BOOST_AUTO_TEST_CASE( single_variant_reader_missing_data )
{
  single_variant_reader_sample_test("testdata/test_variants_missing_data.vcf", vector<string>{"0001A","0002A","0003A","0004A","0005A","0006A","0007A"}, true, 7);       // include 7 samples
  single_variant_reader_sample_test("testdata/test_variants_missing_data.vcf", vector<string>{"0007B", "0008A", "0009A", "0009B"}, false, 7);   // exclude 4 samples

  for (const auto& record : SingleVariantReader{"testdata/test_variants_missing_data.vcf"}){ // include all 11 samples
    BOOST_CHECK_EQUAL(record.n_samples(), 11); 
    const auto gt_for_all_samples = record.genotypes();
    for (const auto& gt_for_single_sample: gt_for_all_samples) {
      BOOST_CHECK(missing(gt_for_single_sample));
      BOOST_CHECK_EQUAL(gt_for_single_sample.size(), 2);
      BOOST_CHECK(missing(gt_for_single_sample.allele_string(0)));
      BOOST_CHECK(missing(gt_for_single_sample.allele_string(1)));
      BOOST_CHECK(missing(gt_for_single_sample.allele_key(0)));
      BOOST_CHECK(missing(gt_for_single_sample.allele_key(1)));
      BOOST_CHECK(missing(gt_for_single_sample[0]));
      BOOST_CHECK(missing(gt_for_single_sample[1]));
    }
  }
}

BOOST_AUTO_TEST_CASE( single_variant_reader_vector )
{
  for (const auto& record : SingleVariantReader{vector<string>{"testdata/test_variants.vcf"}})
    BOOST_CHECK_EQUAL(record.n_samples(), 3);
}

BOOST_AUTO_TEST_CASE( single_variant_reader_vector_all_samples )
{
  for (const auto& record : SingleVariantReader{vector<string>{"testdata/test_variants.vcf"}, vector<string>{}, false})  // include all samples by setting include == false and passing an empty list
    BOOST_CHECK_EQUAL(record.n_samples(), 3);
}

BOOST_AUTO_TEST_CASE( single_variant_reader_vector_too_large )
{
  BOOST_CHECK_THROW((SingleVariantReader{vector<string>{"testdata/test_variants.vcf", "testdata/test_variants.vcf"}}), std::runtime_error);
}

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

BOOST_AUTO_TEST_CASE( multiple_variant_reader_test ) {
  auto truth_index = 0u;
  const auto reader = MultipleVariantReader<MultipleVariantIterator>{{"testdata/test_variants.vcf", "testdata/test_variants.bcf"}, false};
  for (const auto& vec : reader) {
    for (const auto& record : vec) {
      check_variant_basic_api(record, truth_index);
      check_quals_api(record, truth_index);
      check_alt_api(record, truth_index);
      check_filters_api(record, truth_index);
      check_genotype_quals_api(record,truth_index);
      check_phred_likelihoods_api(record, truth_index);
      check_individual_field_api(record, truth_index);
      check_shared_field_api(record, truth_index);
      check_genotype_api(record, truth_index);
    }
    ++truth_index;
  }
}

void multiple_variant_reader_sample_test(const vector<string> samples, const bool include, const int desired_samples) {
  auto filenames = vector<string>{"testdata/test_variants.vcf", "testdata/test_variants.bcf"};

  for (const auto& vec : MultipleVariantReader<MultipleVariantIterator>{filenames, false, samples, include})
    for (const auto& record : vec)
      BOOST_CHECK_EQUAL(record.n_samples(), desired_samples);
}

/*      TODO: Issue #209
BOOST_AUTO_TEST_CASE( multiple_variant_reader_sites_only )
{
  multiple_variant_reader_sample_test(vector<string>{}, true, 0); // exclude all samples (sites-only)
}
*/

BOOST_AUTO_TEST_CASE( multiple_variant_reader_include_all_samples )
{
  multiple_variant_reader_sample_test(vector<string>{}, false, 3);  // include all samples by setting include == false and passing an empty list
}

BOOST_AUTO_TEST_CASE( multiple_variant_reader_including )
{
  multiple_variant_reader_sample_test(vector<string>{"NA12878"}, true, 1); // include only NA12878
  multiple_variant_reader_sample_test(vector<string>{"NA12878", "NA12892"}, true, 2); // include both these samples
}

BOOST_AUTO_TEST_CASE( multiple_variant_reader_excluding )
{
  multiple_variant_reader_sample_test(vector<string>{"NA12891"}, false, 2);  // exclude only NA12891
  multiple_variant_reader_sample_test(vector<string>{"NA12891", "NA12878"}, false, 1);  // exclude both these samples
}

const auto gvcf_truth_ref               = vector<string>{"T", "C", "GG"};
const auto gvcf_truth_chromosome        = vector<uint32_t>{0, 0, 1};
const auto gvcf_truth_alignment_starts  = vector<uint32_t>{10000000, 20000000, 10001000};
const auto gvcf_truth_alignment_stops   = vector<uint32_t>{10000000, 20000123, 10001001};
const auto gvcf_truth_n_alleles         = vector<uint32_t>{2, 2, 2};
const auto gvcf_truth_id                = vector<string>{"db2342", ".", "rs837472"};

BOOST_AUTO_TEST_CASE( gvcf_test ) {
  auto truth_index = 0u;
  const auto reader = SingleVariantReader{"testdata/test.g.vcf"};
  for (const auto& record : reader) {
    BOOST_CHECK_EQUAL(record.ref(), gvcf_truth_ref[truth_index]);
    BOOST_CHECK_EQUAL(record.chromosome(), gvcf_truth_chromosome[truth_index]);
    BOOST_CHECK_EQUAL(record.alignment_start(), gvcf_truth_alignment_starts[truth_index]);
    BOOST_CHECK_EQUAL(record.alignment_stop(), gvcf_truth_alignment_stops[truth_index]);
    BOOST_CHECK_EQUAL(record.n_alleles(), gvcf_truth_n_alleles[truth_index]);
    BOOST_CHECK_EQUAL(record.n_samples(), 3);
    BOOST_CHECK_EQUAL(record.id(), gvcf_truth_id[truth_index]);
    ++truth_index;
  }
}

// TODO?  update to work with VCF GZ
// const auto input_files = vector<string>{"testdata/var_idx/test_variants.bcf", "testdata/var_idx/test_variants_csi.vcf.gz", "testdata/var_idx/test_variants_tabix.vcf.gz"};
const auto indexed_variant_input_files = vector<string>{"testdata/var_idx/test_variants.bcf"};

const auto indexed_variant_chrom_full = vector<string> {"1", "20", "22"};
const auto indexed_variant_bp_full = vector<string> {"1:10000000-10000000", "20:10001000-10001000", "20:10002000-10002000", "20:10003000-10003000", "22:10004000-10004000"};
const auto indexed_variant_chrom_partial = vector<string> {"1"};
const auto indexed_variant_bp_partial = vector<string> {"20:10001000-10001000"};

// empty interval lists and full interval lists
BOOST_AUTO_TEST_CASE( indexed_variant_reader_full_test ) {
  for (const auto filename : indexed_variant_input_files) {
    for (const auto intervals : {vector<string>{}, IndexedVariantIterator::all_intervals, indexed_variant_chrom_full, indexed_variant_bp_full}) {
      auto truth_index = 0u;
      const auto reader = IndexedVariantReader<IndexedVariantIterator>{filename, intervals};
      for (const auto& record : reader) {
        check_variant_basic_api(record, truth_index);
        check_quals_api(record, truth_index);
        check_alt_api(record, truth_index);
        check_filters_api(record, truth_index);
        check_genotype_quals_api(record,truth_index);
        check_phred_likelihoods_api(record, truth_index);
        check_individual_field_api(record, truth_index);
        check_shared_field_api(record, truth_index);
        check_genotype_api(record, truth_index);
        ++truth_index;
      }
      BOOST_CHECK_EQUAL(truth_index, 5u);
    }
  }
}

// one (different) record each
BOOST_AUTO_TEST_CASE( indexed_variant_reader_partial_test ) {
  for (const auto filename : indexed_variant_input_files) {

    const auto intervals1 = indexed_variant_chrom_partial;
    auto truth_index = 0u;
    const auto reader1 = IndexedVariantReader<IndexedVariantIterator>{filename, intervals1};
    for (const auto& record : reader1) {
      BOOST_CHECK_EQUAL(record.ref(), "T");
      BOOST_CHECK_EQUAL(record.chromosome(), 0);
      BOOST_CHECK_EQUAL(record.alignment_start(), 10000000);
      BOOST_CHECK_EQUAL(record.alignment_stop(), 10000000);
      BOOST_CHECK_EQUAL(record.n_alleles(), 2);
      BOOST_CHECK_EQUAL(record.n_samples(), 3);
      BOOST_CHECK_EQUAL(record.id(), "db2342");
      ++truth_index;
    }
    BOOST_CHECK_EQUAL(truth_index, 1u);

    const auto intervals2 = indexed_variant_bp_partial;
    truth_index = 0u;
    const auto reader2 = IndexedVariantReader<IndexedVariantIterator>{filename, intervals2};
    for (const auto& record : reader2) {
      BOOST_CHECK_EQUAL(record.ref(), "GG");
      BOOST_CHECK_EQUAL(record.chromosome(), 1);
      BOOST_CHECK_EQUAL(record.alignment_start(), 10001000);
      BOOST_CHECK_EQUAL(record.alignment_stop(), 10001001);
      BOOST_CHECK_EQUAL(record.n_alleles(), 2);
      BOOST_CHECK_EQUAL(record.n_samples(), 3);
      BOOST_CHECK_EQUAL(record.id(), "rs837472");
      ++truth_index;
    }
    BOOST_CHECK_EQUAL(truth_index, 1u);
  }
}

BOOST_AUTO_TEST_CASE( indexed_variant_reader_move_test ) {
  for (const auto filename : indexed_variant_input_files) {
    auto reader1 = IndexedVariantReader<IndexedVariantIterator>{filename, indexed_variant_chrom_full};
    auto iter1 = reader1.begin();

    // move construct
    auto reader2 = std::move(reader1);
    auto iter2 = reader2.begin();

    // move assign
    reader1 = std::move(reader2);
    auto iter3 = reader1.begin();

    auto rec1 = *iter1;
    auto rec2 = *iter2;
    auto rec3 = *iter3;

    BOOST_CHECK_EQUAL(rec1.alignment_start(), rec2.alignment_start());
    BOOST_CHECK_EQUAL(rec1.alignment_start(), rec3.alignment_start());
  }
 }

