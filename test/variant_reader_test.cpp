#include "variant/variant.h"
#include "variant/variant_reader.h"
#include "variant/multiple_variant_reader.h"
#include "variant/multiple_variant_iterator.h"
#include "variant/indexed_variant_reader.h"
#include "variant/indexed_variant_iterator.h"
#include "variant/synced_variant_reader.h"
#include "variant/synced_variant_iterator.h"
#include "variant/variant_header_builder.h"
#include "exceptions.h"
#include "missing.h"
#include "test_utils.h"

#include <boost/test/unit_test.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <stdexcept>
#include <unordered_set>

using namespace std;
using namespace gamgee;

constexpr auto FLOAT_COMPARISON_THRESHOLD = 0.0001f;

float g_bcf_float_missing;
vector<uint32_t> truth_chromosome;
vector<uint32_t> truth_alignment_starts;
vector<uint32_t> truth_alignment_stops;
vector<uint32_t> truth_n_alleles;
vector<string> truth_filter_name;
vector<uint32_t> truth_filter_size;
vector<float> truth_quals;
vector<string> truth_ref;
vector<vector<string>> truth_alt;
boost::dynamic_bitset<> truth_high_quality_hets;
vector<string> truth_id;
vector<vector<float>> truth_shared_af;
vector<vector<int32_t>> truth_shared_an;
vector<vector<string>> truth_shared_desc;
vector<bool> truth_shared_validated;
vector<int32_t> truth_gq_bt_type;
vector<vector<int32_t>> truth_gq;
vector<vector<vector<float>>> truth_af;
vector<int32_t> truth_pl_bt_type;
vector<vector<vector<int32_t>>> truth_pl;
vector<vector<string>> truth_as;
vector<vector<int32_t>> truth_shared_vlint;
vector<vector<float>> truth_shared_vlfloat;
vector<int32_t> truth_vlint_bt_type;
vector<vector<vector<int32_t>>> truth_vlint;
vector<vector<vector<float>>> truth_vlfloat;

// This gloabl fixture will bear the initialization of all the global variables in this test. 
// The reason we need this is because some of these variables need to be initialized using function
// calls (e.g. bcf_float_set_missing) and having them all initialized as global consts didn't allow me 
// to do this. 
//
// I DO NOT RECOMMEND OR ENDORSE THIS PATTERN
//
// but given that this test file has gone out of control with 
// over 1000 lines that depend on this global initialization pattern, this was the simplest way to dig 
// us out of this hole. In the future, please let's make more use of test facilities (like fixtures) in
// order to avoid getting us deep in the hole again.
struct TruthInitializationFixture {
  TruthInitializationFixture() {
    bcf_float_set_missing(g_bcf_float_missing);
    truth_chromosome        = vector<uint32_t>{0, 1, 1, 1, 2, 2, 2};
    truth_alignment_starts  = vector<uint32_t>{10000000, 10001000, 10002000, 10003000, 10004000, 10005000, 10006000 };
    truth_alignment_stops   = vector<uint32_t>{10000000, 10001001, 10002006, 10003000, 10004002, 10005002, 10006002 };
    truth_n_alleles         = vector<uint32_t>{2, 2, 2, 2, 3, 3, 3};
    truth_filter_name       = vector<string>{"PASS", "PASS", "LOW_QUAL", "NOT_DEFINED", "PASS", "PASS", "PASS"};
    truth_filter_size       = vector<uint32_t>{1,1,1,1,2,1,1};
    truth_quals             = vector<float>{80,8.4,-1,-1,-1,-1,-1};
    truth_ref               = vector<string>{"T", "GG", "TAGTGQA", "A", "GAT","GAT","GAT"};
    truth_alt               = vector< vector<string>> {  { "C" } , {"AA"},  {"T"},  {"AGCT"},  {"G","GATAT"},  {"G","GATAT"}, {"G", "GATAT"}};
    truth_high_quality_hets = boost::dynamic_bitset<>{std::string{"001"}};
    truth_id                = vector<string>{"db2342", "rs837472", ".", ".", ".",".","."};
    truth_shared_af         = vector<vector<float>>{{0.5}, {0.5}, {0.5}, {0.5}, {0.5, 0.0}, {0.5, g_bcf_float_missing}, {0.5, 0.0}};
    truth_shared_an         = vector<vector<int32_t>>{{6}, {6}, {6}, {6}, {6}, {6}, {6}};
    truth_shared_desc       = vector<vector<string>>{{"Test1,Test2"}, {}, {}, {}, {}, {}, {}};
    truth_shared_validated  = vector<bool>{true, false, true, false, false, false, false};
    truth_gq_bt_type        =  vector<int32_t>{BCF_BT_INT16, BCF_BT_INT8, BCF_BT_INT8, BCF_BT_INT8, BCF_BT_INT8, BCF_BT_INT8, BCF_BT_INT8};
    truth_gq                = vector<vector<int32_t>>{{25,12,650}, {35,35,35}, {35,35,35}, {35,35,35}, {35,35,35}, {35,bcf_int32_missing,35}, {35,35,35}};
    truth_af                = vector<vector<vector<float>>> {
      { {3.1,2.2}, 		   {3.1,2.2}, {3.1,2.2} }, 
        { {3.1,2.2}, 		   {3.1,2.2}, {3.1,2.2} }, 
        { {3.1,2.2}, 		   {3.1,2.2}, {3.1,2.2} }, 
        { {3.1,2.2}, 		   {3.1,2.2}, {3.1,2.2} }, 
        { {3.1,2.2}, 		   {3.1,2.2}, {3.1,2.2} }, 
        { {3.1,g_bcf_float_missing}, {3.1,2.2}, {3.1,2.2} },
        { {3.1,2.2}, 		   {3.1,2.2}, {3.1,2.2} }
    };
    truth_pl_bt_type = vector<int32_t>{ BCF_BT_INT8, BCF_BT_INT8, BCF_BT_INT32, BCF_BT_INT8, BCF_BT_INT8, BCF_BT_INT8, BCF_BT_INT8 };
    truth_pl                = vector<vector<vector<int32_t>>>{
      {{10,0,100      }, {0,10,1000      }, {10,100,0}      },
        {{10,0,100      }, {0,10,100       }, {10,100,0}      },
        {{10,0,100      }, {0,10,2000000000}, {10,100,0}      },
        {{10,0,100      }, {0,10,100       }, {10,100,0}      },
        {{10,0,100,2,4,8}, {0,10,100,2,4,8 }, {10,100,0,2,4,8}},
        {{10,0,100,bcf_int32_missing,4,bcf_int32_missing}, {0,10,100,2,4,8 }, {10,100,0,2,4,8}},
        {{10,0,100,2,4,8}, {0,10,100,2,4,8 }, {10,100,0,2,4,8}}
    };
    truth_as                = vector<vector<string>>{ 
      {"ABA","CA","XISPAFORRO"}, 
        {"ABA","ABA","ABA"}, 
        {"ABA","ABA","ABA"}, 
        {"ABA","ABA","ABA"}, 
        {"ABA","ABA","ABA"}, 
        {"ABA","ABA","."},
        {"ABA","ABA","ABA"}
    };
    truth_shared_vlint = vector<vector<int32_t>>{{}, {}, {}, {}, {}, {}, {27, 57, 122}};
    truth_shared_vlfloat = vector<vector<float>>{{}, {}, {}, {}, {}, {}, {5.3, -127.65, 24245.9, 435.78}};
    truth_vlint_bt_type = vector<int32_t>{ -1, -1, -1, -1, -1, -1, BCF_BT_INT16 };
    truth_vlint                = vector<vector<vector<int32_t>>>{
      {{ }, { }, { }},
	{{ }, { }, { }},
	{{ }, { }, { }},
	{{ }, { }, { }},
	{{ }, { }, { }},
	{{ }, { }, { }},
	{{4,56,21}, {36,bcf_int32_missing,1024,5196}, {2,1}}
    };
    truth_vlfloat                = vector<vector<vector<float>>>{
      {{ }, { }, { }},
	{{ }, { }, { }},
	{{ }, { }, { }},
	{{ }, { }, { }},
	{{ }, { }, { }},
	{{ }, { }, { }},
	{{4.5}, {9.25,-15.125}, {100.5,g_bcf_float_missing,-92.75,-16345.25}}
    };
  }
};

BOOST_GLOBAL_FIXTURE(TruthInitializationFixture);

/*
 * @brief given the data type (BCF_BT_*), returns the int32_t with the value of bcf_*_vector_end
 * @param bcf_bt_type type (BCF_BT_*) of the data
 * @return bcf_*_vector_end
 */
int32_t bcf_vector_end_val(int32_t bcf_bt_type)
{
  auto vector_end_val = 0;
  switch(bcf_bt_type)
  {
    case BCF_BT_INT8:
      vector_end_val = bcf_int8_vector_end;
      break;
    case BCF_BT_INT16:
      vector_end_val = bcf_int16_vector_end;
      break;
    case BCF_BT_INT32:
      vector_end_val = bcf_int32_vector_end;
      break;
    case BCF_BT_FLOAT:
      vector_end_val = bcf_float_vector_end;
      break;
    default:
      BOOST_CHECK_MESSAGE(0, "Unknown VCF integer type "<< bcf_bt_type);
      return -1;
  }
  return vector_end_val;
}

/*
 * @brief given the int data type (BCF_BT_INT*) and int32 representation of the datatype, returns the string corresponding to 
 * the data type. Note that functions like shared_field_as_integer return int32 representation of data, even though the underlying data
 * type might be int16/int8 etc. Hence, shared_field_as_integer would convert bcf_int8_missing to bcf_int32_missing so that
 * missing functions in the rest of gamgee work correctly. However, a string comparison check in the tests here would fail 
 * since to_string(bcf_int8_missing) != to_string(bcf_int32_missing). This function re-converts 'special' values back to their
 * original data types and then converts to string.
 * @note the function checks whether it is a special value - anything <= bcf_int32_vector_end is a special value- if yes, compute
 * the real special value corresponding to the data type bcf_bt_type. The 'real' special value is computed by obtaining the 
 * bcf_*_vector_end value and subtracting from it 'diff', where diff is the difference between v and the int32 special value
 * bcf_int32_vector_end
 * @param v int32_t value obtained from *_as_integer gamgee function
 * @param bcf_bt_type int type (BCF_BT_INT*) of the data
 * @return string representation of value
 */
string bcf_int32_to_string(int32_t v, int32_t bcf_bt_type)
{
  if(v <= bcf_int32_vector_end)
  {
    auto vector_end_val = bcf_vector_end_val(bcf_bt_type);    
    auto diff = bcf_int32_vector_end - v; //either 0 or 1
    return to_string(vector_end_val - diff);
  }
  else
    return to_string(v);
}
/*
 * @brief given the int data type (BCF_BT_INT*) and int32 representation of the datatype, returns the float corresponding to 
 * the data type. Note that functions like shared_field_as_integer return int32 representation of data, even though the underlying data
 * type might be int16/int8 etc. Hence, shared_field_as_integer would convert bcf_int8_missing to bcf_int32_missing so that
 * missing functions in the rest of gamgee work correctly. However, a float comparison check in the tests here would fail 
 * since float(bcf_int8_missing) != float(bcf_int32_missing). This function re-converts 'special' values back to their
 * original data types and then converts to float.
 * @note the function checks whether it is a special value - anything <= bcf_int32_vector_end is a special value- if yes, compute
 * the real special value corresponding to the data type bcf_bt_type. The 'real' special value is computed by obtaining the 
 * bcf_*_vector_end value and subtracting from it 'diff', where diff is the difference between v and the int32 special value
 * bcf_int32_vector_end
 * @param v int32_t value obtained from *_as_integer gamgee function
 * @param bcf_bt_type type (BCF_BT_*) of the data
 * @return string representation of value
 */
float bcf_int32_to_float(int32_t v, int32_t bcf_bt_type)
{
  if(v <= bcf_int32_vector_end)
  {
    auto diff = bcf_int32_vector_end - v; //either 0 or 1
    auto pack_value = bcf_float_vector_end - diff;
    auto tmp = 0.0f;
    bcf_float_set(&tmp, pack_value);
    return tmp;
  }
  else
    return float(v);
}
/*
 * @brief floating point comparison for VCF fields
 * v1, v2 can be NaNs since the bcf_float_missing and vector_end values are NaNs. A simple comparison
 * of v1 and v2 would fail since (nan == nan) is always false. This is fixed by explicitly comparing both values to bcf_float_missing and
 * vector_end i.e. v1 and v2 are equal if ANY ONE of the following conditions is true
 * (a) both are bcf_float_missing OR
 * (b) both are bcf_float_vector_end OR
 * (c) |v1-v2|<threshold
 */
bool bcf_compare_float(const float v1, const float v2, const float threshold)
{
  auto abs_diff = fabsf(v1-v2);
  auto rel_diff = (v1 != 0) ? abs_diff/fabsf(v1) : 0;
  return ((bcf_float_is_missing(v1) && bcf_float_is_missing(v2)) || (bcf_float_is_vector_end(v1) && bcf_float_is_vector_end(v2))
	|| (abs_diff <= threshold || rel_diff <= threshold));
}


boost::dynamic_bitset<> high_qual_hets(const Variant& record) {  // filter all hets that have GQ > 20
  const auto genotypes = record.genotypes(); // a "vector-like" with the genotypes of all samples in this record
  const auto gqs = record.integer_individual_field("GQ"); // a "vector-like" with all the GQs of all samples in this record
  const auto hets = Variant::select_if(genotypes.begin(), genotypes.end(), [](const Genotype& g) { return g.het(); }); // returns a bit set with all hets marked with 1's
  const auto pass_gqs = Variant::select_if(gqs.begin(), gqs.end(), [](const IndividualFieldValue<int32_t>& gq) { return gq[0] > 20; }); // returns a bit set with every sample with gq > 20 marked with 1's
  return hets & pass_gqs; // returns a bit set with all the samples that are het and have gq > 20
}

void check_variant_basic_api(const Variant& record, const uint32_t truth_index) {
  BOOST_CHECK_EQUAL(record.ref(), truth_ref[truth_index]);
  BOOST_CHECK_EQUAL(record.chromosome(), truth_chromosome[truth_index]);
  BOOST_CHECK_EQUAL(record.alignment_start(), truth_alignment_starts[truth_index]);
  BOOST_CHECK_EQUAL(record.alignment_stop(), truth_alignment_stops[truth_index]);
  BOOST_CHECK_EQUAL(record.n_alleles(), truth_n_alleles[truth_index]);
  BOOST_CHECK_EQUAL(record.n_samples(), 3u);
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
  BOOST_CHECK_EQUAL(count_if(filters.begin(), filters.end(), [](const string& x){return x == x;}), truth_filter_size[truth_index]); // checking iteration in functional style
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
    BOOST_CHECK_EQUAL(1, count_if(sample_pl.begin(), sample_pl.end(), [](const int x) { return x == 0; })); // testing std library functional call at the sample value level
    for(const auto& value_pl : sample_pl)  // testing for-each iteration at the sample value level
      BOOST_CHECK_EQUAL(value_pl, value_pl); // I just want to make sure that this level of iteration is supported, the values don't matter anymore at this point
  }
}      

void check_variable_length_field_api(const Variant& record, const uint32_t truth_index) {
  const auto header = record.header();
  const auto vlint_int = record.individual_field_as_integer("VLINT");
  const auto vlint_float = record.individual_field_as_float("VLINT");
  const auto vlint_string = record.individual_field_as_string("VLINT");
  const auto vlfloat_int = record.individual_field_as_integer("VLFLOAT");
  const auto vlfloat_float = record.individual_field_as_float("VLFLOAT");
  const auto vlfloat_string = record.individual_field_as_string("VLFLOAT");
  const auto vlint_idx = header.field_index("VLINT");
  const auto vlfloat_idx = header.field_index("VLFLOAT");
  const auto vlint_int_idx    = record.individual_field_as_integer(vlint_idx);
  const auto vlint_float_idx  = record.individual_field_as_float(vlint_idx);
  const auto vlint_string_idx = record.individual_field_as_string(vlint_idx); 
  const auto vlfloat_int_idx    = record.individual_field_as_integer(vlfloat_idx);
  const auto vlfloat_float_idx  = record.individual_field_as_float(vlfloat_idx);
  const auto vlfloat_string_idx = record.individual_field_as_string(vlfloat_idx); 
  BOOST_CHECK(vlint_int == record.integer_individual_field("VLINT"));
  BOOST_CHECK(vlfloat_float == record.float_individual_field("VLFLOAT"));
  if(!vlfloat_float.empty())
    BOOST_CHECK(vlfloat_float != record.individual_field_as_float("VLINT"));
  for(auto i=0u; i != record.n_samples(); ++i) {
    if(!vlint_int.empty())
    {
      auto j = 0u;
      //the vector shoud have at least as many values as the truth
      BOOST_CHECK(truth_vlint[truth_index][i].size() <= vlint_int[i].size());
      for (auto vlint_value : vlint_int[i])
      {
	BOOST_CHECK_MESSAGE(!utils::bcf_is_vector_end_value<int32_t>(vlint_value), "Saw vector end value : hex value "<<hex<<vlint_value);
	BOOST_CHECK(j < truth_vlint[truth_index][i].size());
	if(j >= truth_vlint[truth_index][i].size())
	  break;
	BOOST_CHECK_EQUAL(vlint_value, truth_vlint[truth_index][i][j]);
	BOOST_CHECK_EQUAL(vlint_int_idx[i][j], truth_vlint[truth_index][i][j]);
	BOOST_CHECK(bcf_compare_float(vlint_float[i][j], bcf_int32_to_float(truth_vlint[truth_index][i][j], truth_vlint_bt_type[truth_index]), FLOAT_COMPARISON_THRESHOLD));
	BOOST_CHECK(bcf_compare_float(vlint_float_idx[i][j], bcf_int32_to_float(truth_vlint[truth_index][i][j], truth_vlint_bt_type[truth_index]), FLOAT_COMPARISON_THRESHOLD));
	BOOST_CHECK_EQUAL(vlint_string[i][j], bcf_int32_to_string(truth_vlint[truth_index][i][j], truth_vlint_bt_type[truth_index]));
	BOOST_CHECK_EQUAL(vlint_string_idx[i][j], bcf_int32_to_string(truth_vlint[truth_index][i][j], truth_vlint_bt_type[truth_index]));
	++j;
      }
    }
    if(!vlfloat_float.empty())
    {
      auto j = 0u;
      //the vector shoud have at least as many values as the truth
      BOOST_CHECK(truth_vlfloat[truth_index][i].size() <= vlfloat_float[i].size());
      for (auto vlfloat_value : vlfloat_float[i])
      {
	BOOST_CHECK_MESSAGE(!utils::bcf_is_vector_end_value<float>(vlfloat_value), "Saw vector end value "<<vlfloat_value);
	BOOST_CHECK(j < truth_vlfloat[truth_index][i].size());
	if(j >= truth_vlfloat[truth_index][i].size())
	  break;
	BOOST_CHECK_EQUAL(vlfloat_int[i][j], int32_t(truth_vlfloat[truth_index][i][j]));
	BOOST_CHECK_EQUAL(vlfloat_int_idx[i][j], int32_t(truth_vlfloat[truth_index][i][j]));
	BOOST_CHECK(bcf_compare_float(vlfloat_float[i][j], truth_vlfloat[truth_index][i][j], FLOAT_COMPARISON_THRESHOLD));
	BOOST_CHECK(bcf_compare_float(vlfloat_float_idx[i][j], truth_vlfloat[truth_index][i][j], FLOAT_COMPARISON_THRESHOLD));
	//For large floating point values, the to_string does not always match
	//For example, 16457.24 and 16457.09 are relatively close but not absolutely close wrt threshold
	//BOOST_CHECK_EQUAL(vlfloat_string[i][j], to_string(truth_vlfloat[truth_index][i][j]));
	//BOOST_CHECK_EQUAL(vlfloat_string_idx[i][j], to_string(truth_vlfloat[truth_index][i][j]));
	++j;
      }
    }
  }
}

//check operator[] for shared and individual iterators
void check_operator_index_for_iterators(const Variant& record)
{
  //Check operator[] for SharedFieldIterator
  const auto af_float  = record.individual_field_as_float("AF");
  auto af_float_iter = af_float.begin();
  for(auto af_float_element : af_float)
  {
    //test offset 0 for iter
    BOOST_CHECK(af_float_element == af_float_iter[0]);
    ++af_float_iter;
    //test non-0 offset for iter
    if(af_float.size() > 0)
    {
      auto af_float_iter_1 = af_float.begin();
      ++af_float_iter_1;
      auto j = 0u;
      for(auto af_float_element : af_float)
      {
	if(j > 0)
	  BOOST_CHECK(af_float_element == af_float_iter_1[j-1]);
	++j;
      }
    }
  }
  //Check operator[] for IndividualFieldIterator & IndividualFieldValueIterator
  const auto pl_int    = record.individual_field_as_integer("PL");
  auto pl_int_iter = pl_int.begin();
  //test offset 0 for iter
  for(auto pl_int_sample : pl_int)
  {
    BOOST_CHECK(pl_int_sample == pl_int_iter[0]);
    ++pl_int_iter;
    //test offset 0 for iter
    auto pl_int_sample_iter = pl_int_sample.begin();
    for(auto pl_int_sample_element : pl_int_sample)
    {
      BOOST_CHECK(pl_int_sample_element == pl_int_sample_iter[0]);
      ++pl_int_sample_iter;
    }
    //test non-0 offset for iter
    if(pl_int_sample.size() > 0)
    {
      auto pl_int_sample_iter_1 = pl_int_sample.begin();
      ++pl_int_sample_iter_1;
      auto j = 0u;
      for(auto pl_int_sample_element : pl_int_sample)
      {
	if(j > 0)
	  BOOST_CHECK(pl_int_sample_element == pl_int_sample_iter_1[j-1]);
	++j;
      }
    }
  }
  //test non-0 offset for iter
  if(pl_int.size() > 0)
  {
    auto pl_int_iter_1 = pl_int.begin();
    ++pl_int_iter_1;
    auto j = 0u;
    for(auto pl_int_sample : pl_int)
    {
      if(j > 0)
	BOOST_CHECK(pl_int_sample == pl_int_iter_1[j-1]);
      ++j;
    }
  }
}

void check_field_iterator_inc_dec(const Variant& record)
{
  const auto pl_int  = record.individual_field_as_integer("PL");
  for(auto pl : pl_int)
  {
    for(auto pl_iter = pl.begin();pl_iter != pl.end();)
    {
      auto curr_elem = *pl_iter;
      BOOST_CHECK_EQUAL(*(pl_iter++), curr_elem);
    }
    auto pl_next_iter=pl.begin();
    ++pl_next_iter;
    for(auto pl_iter = pl.begin();pl_next_iter != pl.end();++pl_next_iter)
    {
      auto pl_tmp_iter = pl_next_iter;
      BOOST_CHECK_EQUAL(*(pl_tmp_iter--), *pl_next_iter);
      BOOST_CHECK_EQUAL(*pl_tmp_iter, *pl_iter);
      ++pl_tmp_iter;
      BOOST_CHECK_EQUAL(*(--pl_tmp_iter), *pl_iter);
      BOOST_CHECK_EQUAL(*(++pl_iter), *pl_next_iter);
    }
  }

  const auto vlint  = record.shared_field_as_integer("VLINT");
  for(auto vlint_iter = vlint.begin();vlint_iter != vlint.end();)
  {
    auto curr_elem = *vlint_iter;
    BOOST_CHECK_EQUAL(*(vlint_iter++), curr_elem);
  }
  auto vlint_next_iter=vlint.begin();
  ++vlint_next_iter;
  for(auto vlint_iter = vlint.begin();vlint_next_iter != vlint.end();++vlint_next_iter)
  {
    auto vlint_tmp_iter = vlint_next_iter;
    BOOST_CHECK_EQUAL(*(vlint_tmp_iter--), *vlint_next_iter);
    BOOST_CHECK_EQUAL(*vlint_tmp_iter, *vlint_iter);
    ++vlint_tmp_iter;
    BOOST_CHECK_EQUAL(*(--vlint_tmp_iter), *vlint_iter);
    BOOST_CHECK_EQUAL(*(++vlint_iter), *vlint_next_iter);
  }
}

void check_out_of_bound_exceptions(const Variant& record) {
  const auto vlint_shared = record.integer_shared_field("VLINT");
  const auto af_float_shared  = record.shared_field_as_float("AF");
  if(!vlint_shared.empty())
    BOOST_CHECK_THROW(vlint_shared[vlint_shared.size()], out_of_range);
  if(!af_float_shared.empty())
    BOOST_CHECK_THROW(af_float_shared[record.n_alleles()], out_of_range);

  const auto as_string = record.individual_field_as_string("AS");
  const auto pl_int    = record.individual_field_as_integer("PL");
  const auto vlfloat  = record.individual_field_as_float("VLFLOAT");
  auto n_alleles = record.n_alleles();
  auto num_gts = (n_alleles*(n_alleles+1))/2;
  auto n_samples = record.n_samples(); 
  for(auto i=0u; i != record.n_samples(); ++i) {
    if(as_string.empty())
    {
      BOOST_CHECK_THROW(as_string[i][0], out_of_range);
    }
    else
    {
      BOOST_CHECK_THROW(as_string[n_samples][0], out_of_range);
      BOOST_CHECK_THROW(as_string[i][1], out_of_range);
    }
    if(pl_int.empty())
    {
      BOOST_CHECK_THROW(pl_int[i][0], out_of_range);
    }
    else
    {
      BOOST_CHECK_THROW(pl_int[n_samples][0], out_of_range);
      BOOST_CHECK_THROW(pl_int[i][num_gts], out_of_range);
    }
    if(vlfloat.empty())
    {
      BOOST_CHECK_THROW(vlfloat[i][0], out_of_range);
    }
    else
    {
      BOOST_CHECK_THROW(vlfloat[n_samples][0], out_of_range);
      BOOST_CHECK_THROW(vlfloat[i][vlfloat[i].size()], out_of_range);
    }
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

  // check types
  BOOST_CHECK_EQUAL(record.header().individual_field_type("GQ"), BCF_HT_INT);
  BOOST_CHECK_EQUAL(record.header().individual_field_type("AF"), BCF_HT_REAL);
  BOOST_CHECK_EQUAL(record.header().individual_field_type("AS"), BCF_HT_STR);
  BOOST_CHECK_EQUAL(record.header().individual_field_type(header.field_index("GQ")), BCF_HT_INT);
  BOOST_CHECK_EQUAL(record.header().individual_field_type(header.field_index("AF")), BCF_HT_REAL);
  BOOST_CHECK_EQUAL(record.header().individual_field_type(header.field_index("AS")), BCF_HT_STR);

  for(auto i=0u; i != record.n_samples(); ++i) {
    BOOST_CHECK_EQUAL(gq_int[i][0], truth_gq[truth_index][i]);
    BOOST_CHECK_EQUAL(gq_int_idx[i][0], truth_gq[truth_index][i]);
    BOOST_CHECK_MESSAGE(bcf_compare_float(gq_float[i][0],
	  bcf_int32_to_float(truth_gq[truth_index][i], truth_gq_bt_type[truth_index]), FLOAT_COMPARISON_THRESHOLD),
	"["<<gq_float[i][0]<<" != "<<truth_gq[truth_index][i]<<"]");
    BOOST_CHECK_MESSAGE(bcf_compare_float(gq_float_idx[i][0],
	  bcf_int32_to_float(truth_gq[truth_index][i], truth_gq_bt_type[truth_index]), FLOAT_COMPARISON_THRESHOLD),
	"["<<gq_float_idx[i][0]<<" != "<<truth_gq[truth_index][i]<<"]");
    BOOST_CHECK_EQUAL(gq_string[i][0], bcf_int32_to_string(truth_gq[truth_index][i], truth_gq_bt_type[truth_index]));
    BOOST_CHECK_EQUAL(gq_string_idx[i][0], bcf_int32_to_string(truth_gq[truth_index][i], truth_gq_bt_type[truth_index]));
    BOOST_REQUIRE_EQUAL(af_int[i].size(), truth_af[truth_index][i].size()); // require otherwise next line may segfault
    BOOST_REQUIRE_EQUAL(af_int_idx[i].size(), truth_af[truth_index][i].size()); // require otherwise next line may segfault
    for (auto j=0u; j!= af_int[i].size(); ++j) {
      BOOST_CHECK_EQUAL(af_int[i][j], int32_t(truth_af[truth_index][i][j]));
      BOOST_CHECK_EQUAL(af_int_idx[i][j], int32_t(truth_af[truth_index][i][j]));
      BOOST_CHECK(bcf_compare_float(af_float[i][j], truth_af[truth_index][i][j], FLOAT_COMPARISON_THRESHOLD));
      BOOST_CHECK(bcf_compare_float(af_float_idx[i][j], truth_af[truth_index][i][j], FLOAT_COMPARISON_THRESHOLD));
      BOOST_CHECK_EQUAL(af_string[i][j], to_string(truth_af[truth_index][i][j]));
      BOOST_CHECK_EQUAL(af_string_idx[i][j], to_string(truth_af[truth_index][i][j]));
    }
    BOOST_REQUIRE_EQUAL(pl_int[i].size(), truth_pl[truth_index][i].size()); // require otherwise next line may segfault
    BOOST_REQUIRE_EQUAL(pl_int_idx[i].size(), truth_pl[truth_index][i].size()); // require otherwise next line may segfault
    for (auto j=0u; j!= pl_int[i].size(); ++j) { 
      BOOST_CHECK_EQUAL(pl_int[i][j], truth_pl[truth_index][i][j]);
      BOOST_CHECK_EQUAL(pl_int_idx[i][j], truth_pl[truth_index][i][j]);
      BOOST_CHECK(bcf_compare_float(pl_float[i][j], bcf_int32_to_float(truth_pl[truth_index][i][j], truth_pl_bt_type[truth_index]), FLOAT_COMPARISON_THRESHOLD));
      BOOST_CHECK(bcf_compare_float(pl_float_idx[i][j], bcf_int32_to_float(truth_pl[truth_index][i][j], truth_pl_bt_type[truth_index]), FLOAT_COMPARISON_THRESHOLD));
      BOOST_CHECK_EQUAL(pl_string[i][j], bcf_int32_to_string(truth_pl[truth_index][i][j], truth_pl_bt_type[truth_index]));
      BOOST_CHECK_EQUAL(pl_string_idx[i][j], bcf_int32_to_string(truth_pl[truth_index][i][j], truth_pl_bt_type[truth_index]));
    }
    BOOST_CHECK_EQUAL(as_string[i][0], truth_as[truth_index][i]);
    BOOST_CHECK_EQUAL(as_string_idx[i][0], truth_as[truth_index][i]);
  }
  BOOST_CHECK_EQUAL(gq_int.front()[0], truth_gq[truth_index][0]);
  BOOST_CHECK_EQUAL(gq_int.back()[0], truth_gq[truth_index][gq_int.size() - 1]);
  BOOST_CHECK_EQUAL(gq_int[0].front(), truth_gq[truth_index][0]);
  BOOST_CHECK_EQUAL(gq_int[0].back(), truth_gq[truth_index][0]);
  BOOST_CHECK_EQUAL(pl_int.front()[0], truth_pl[truth_index][0][0]);
  BOOST_CHECK_EQUAL(pl_int.back()[0], truth_pl[truth_index][pl_int.size() - 1][0]);
  BOOST_CHECK_EQUAL(pl_int[0].front(), truth_pl[truth_index][0][0]);
  BOOST_CHECK_EQUAL(pl_int[0].back(), truth_pl[truth_index][0][pl_int[0].size() - 1]);
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
  check_variable_length_field_api(record, truth_index);
  check_operator_index_for_iterators(record);
  check_out_of_bound_exceptions(record);
  check_field_iterator_inc_dec(record);
}

void check_shared_field_api(const Variant& record, const uint32_t truth_index) {
  const auto header = record.header();
  BOOST_CHECK_EQUAL(record.boolean_shared_field("VALIDATED"), truth_shared_validated[truth_index]);
  const auto an = record.integer_shared_field("AN"); // test the string based api
  BOOST_CHECK_EQUAL_COLLECTIONS(an.begin(), an.end(), truth_shared_an[truth_index].begin(), truth_shared_an[truth_index].end());
  const auto an_idx = record.integer_shared_field(header.field_index("AN")); // test the index based api
  BOOST_CHECK_EQUAL_COLLECTIONS(an_idx.begin(), an_idx.end(), truth_shared_an[truth_index].begin(), truth_shared_an[truth_index].end());
  const auto af = record.float_shared_field("AF");
  int af_index = 0;
  for(auto af_value : af)
  {
    BOOST_CHECK_MESSAGE(bcf_compare_float(af_value, truth_shared_af[truth_index][af_index], FLOAT_COMPARISON_THRESHOLD),
	"["<<af_value<<" != "<<truth_shared_af[truth_index][af_index]<<"]");
    ++af_index;
  }
  //More C++ 2011 friendly - but mismatch with const iter makes it impossible to use without larger code changes
  //for_each(boost::make_zip_iterator(boost::make_tuple(af.begin(), truth_shared_af[truth_index].cbegin())),
      //boost::make_zip_iterator(boost::make_tuple(af.end(), truth_shared_af[truth_index].cend())), [](auto tuple)
      //{ bcf_compare_float(tuple.get<0>(),
	//tuple.get<1>(), FLOAT_COMPARISON_THRESHOLD);  } );
  const auto af_idx = record.float_shared_field(header.field_index("AF")); // test the index based api
  af_index = 0;
  for(auto af_value : af_idx)
  {
    BOOST_CHECK_MESSAGE(bcf_compare_float(af_value, truth_shared_af[truth_index][af_index], FLOAT_COMPARISON_THRESHOLD),
	"["<<af_value<<" != "<<truth_shared_af[truth_index][af_index]<<"]");
    ++af_index;
  }
  const auto vlint = record.integer_shared_field("VLINT");
  if(!vlint.empty())
  {
    BOOST_CHECK_EQUAL(vlint.size(), truth_shared_vlint[truth_index].size());
    BOOST_CHECK_EQUAL_COLLECTIONS(vlint.begin(), vlint.end(), truth_shared_vlint[truth_index].begin(), truth_shared_vlint[truth_index].end());
  }
  const auto vlfloat = record.float_shared_field("VLFLOAT");
  if(!vlfloat.empty())
  {
    BOOST_CHECK_EQUAL(vlfloat.size(), truth_shared_vlfloat[truth_index].size());
    BOOST_CHECK_EQUAL_COLLECTIONS(vlfloat.begin(), vlfloat.end(), truth_shared_vlfloat[truth_index].begin(), truth_shared_vlfloat[truth_index].end());
  }
  const auto desc = record.string_shared_field("DESC");
  BOOST_CHECK_EQUAL_COLLECTIONS(desc.begin(), desc.end(), truth_shared_desc[truth_index].begin(), truth_shared_desc[truth_index].end());
  const auto desc_idx = record.string_shared_field(header.field_index("DESC")); // test the index based api
  BOOST_CHECK_EQUAL_COLLECTIONS(desc_idx.begin(), desc_idx.end(), truth_shared_desc[truth_index].begin(), truth_shared_desc[truth_index].end());
  // check front() and back()
  BOOST_CHECK_EQUAL(an.front(), truth_shared_an[truth_index][0]);
  BOOST_CHECK_EQUAL(an.back(), truth_shared_an[truth_index][an.size() - 1]);
  // check non-existing missing values
  BOOST_CHECK(missing(record.boolean_shared_field("NON_EXISTING")));
  BOOST_CHECK(missing(record.integer_shared_field("NON_EXISTING")));
  BOOST_CHECK(missing(record.float_shared_field("NON_EXISTING")));  
  BOOST_CHECK(missing(record.string_shared_field("NON_EXISTING"))); 
  BOOST_CHECK(missing(record.boolean_shared_field(-1)));
  BOOST_CHECK(missing(record.integer_shared_field(-1)));
  BOOST_CHECK(missing(record.float_shared_field(-1)));  
  BOOST_CHECK(missing(record.string_shared_field(-1))); 
  // type checking
  BOOST_CHECK_EQUAL(record.header().shared_field_type("VALIDATED"), BCF_HT_FLAG);
  BOOST_CHECK_EQUAL(record.header().shared_field_type("AN"), BCF_HT_INT);
  BOOST_CHECK_EQUAL(record.header().shared_field_type("AF"), BCF_HT_REAL);
  BOOST_CHECK_EQUAL(record.header().shared_field_type("DESC"), BCF_HT_STR);
  BOOST_CHECK_EQUAL(record.header().shared_field_type(header.field_index("VALIDATED")), BCF_HT_FLAG);
  BOOST_CHECK_EQUAL(record.header().shared_field_type(header.field_index("AN")), BCF_HT_INT);
  BOOST_CHECK_EQUAL(record.header().shared_field_type(header.field_index("AF")), BCF_HT_REAL);
  BOOST_CHECK_EQUAL(record.header().shared_field_type(header.field_index("DESC")), BCF_HT_STR);
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
    BOOST_CHECK_EQUAL(gt_for_single_sample.size(), 2u);
    const auto alleles = gt_for_single_sample.allele_strings();
    const auto allele_keys = gt_for_single_sample.allele_keys();
    BOOST_CHECK_EQUAL(gt_for_single_sample.allele_string(0), alleles[0]);
    BOOST_CHECK_EQUAL(gt_for_single_sample.allele_string(1), alleles[1]);
    BOOST_CHECK_EQUAL(gt_for_single_sample.allele_key(0), allele_keys[0]);
    BOOST_CHECK_EQUAL(gt_for_single_sample.allele_key(1), allele_keys[1]);
    BOOST_CHECK_EQUAL(gt_for_single_sample[0], allele_keys[0]);
    BOOST_CHECK_EQUAL(gt_for_single_sample[1], allele_keys[1]);
    BOOST_CHECK(!missing(gt_for_single_sample));
    BOOST_CHECK_EQUAL(gt_for_single_sample.size(), 2u);
    if (gt_for_single_sample.het()) {
      BOOST_CHECK(!gt_for_single_sample.hom_ref());
      BOOST_CHECK(!gt_for_single_sample.hom_var());
      BOOST_CHECK_NE(allele_keys[0], allele_keys[1]);
      BOOST_CHECK_NE(alleles[0], alleles[1]);
      if (truth_index != 4) {
        BOOST_CHECK(!gt_for_single_sample.non_ref_het());
        BOOST_CHECK_EQUAL(allele_keys[0], 0);
        BOOST_CHECK_EQUAL(allele_keys[0], 0);
        BOOST_CHECK_EQUAL(allele_keys[1], 1);
        BOOST_CHECK_EQUAL(alleles[0], truth_ref[truth_index]);
        BOOST_CHECK_EQUAL(alleles[1], truth_alt[truth_index][0]);
        BOOST_CHECK_EQUAL(gt_for_single_sample.fast_diploid_key_generation(), 0x00000001u);
      } else {
        BOOST_CHECK(gt_for_single_sample.non_ref_het());
        BOOST_CHECK_EQUAL(allele_keys[0], 1);
        BOOST_CHECK_EQUAL(allele_keys[1], 2);
        BOOST_CHECK_NE(alleles[0], alleles[1]);
        BOOST_CHECK_EQUAL(alleles[0], truth_alt[truth_index][0]);
        BOOST_CHECK_EQUAL(alleles[1], truth_alt[truth_index][1]);
        BOOST_CHECK_EQUAL(gt_for_single_sample.fast_diploid_key_generation(), 0x00010002u);
      }
    }
    if (gt_for_single_sample.hom_ref()) {
      BOOST_CHECK(!gt_for_single_sample.het());
      BOOST_CHECK(!gt_for_single_sample.non_ref_het());
      BOOST_CHECK(!gt_for_single_sample.hom_var());
      BOOST_CHECK_EQUAL(allele_keys[0], allele_keys[1]);
      BOOST_CHECK_EQUAL(alleles[0], alleles[1]);
      BOOST_CHECK_EQUAL(allele_keys[0], 0);
      BOOST_CHECK_EQUAL(allele_keys[1], 0);
      BOOST_CHECK_EQUAL(alleles[0], truth_ref[truth_index]);
      BOOST_CHECK_EQUAL(alleles[1], truth_ref[truth_index]);
      BOOST_CHECK_EQUAL(gt_for_single_sample.fast_diploid_key_generation(), 0x00000000u);
    }
    if (gt_for_single_sample.hom_var()) {
      BOOST_CHECK(!gt_for_single_sample.het());
      BOOST_CHECK(!gt_for_single_sample.non_ref_het());
      BOOST_CHECK(!gt_for_single_sample.hom_ref());
      BOOST_CHECK_EQUAL(allele_keys[0], allele_keys[1]);
      BOOST_CHECK_EQUAL(alleles[0], alleles[1]);
      BOOST_CHECK_EQUAL(allele_keys[0], 1);
      BOOST_CHECK_EQUAL(allele_keys[1], 1);
      BOOST_CHECK_EQUAL(alleles[0], truth_alt[truth_index][0]);
      BOOST_CHECK_EQUAL(alleles[1], truth_alt[truth_index][0]);
      BOOST_CHECK_EQUAL(gt_for_single_sample.fast_diploid_key_generation(), 0x00010001u);
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

void check_all_apis(const Variant& record, const uint32_t truth_index) {
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

void generic_variant_reader_test(const std::function<void(const Variant&, const uint32_t)>& fun) {
  for (const auto& filename : {"testdata/test_variants.vcf", "testdata/test_variants.bcf", "testdata/test_variants.vcf.gz"}) {
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
  const auto truth_missing = vector<bool>{false, false, true, true, true, true, true};
  auto i = 0u;
  for (const auto& record : SingleVariantReader{"testdata/test_variants.vcf"})
    BOOST_CHECK_EQUAL(missing(record.id()), truth_missing[i++]); // check that the missing and non-missing values in the vcf actually return the right missingness
  BOOST_CHECK(missing(""));                                      // check that an empty string is missing
  BOOST_CHECK(!missing("arBItrary_vaLue"));                      // check that a non-empty string is not missing
}

void single_variant_reader_sample_test(const string filename, const vector<string> samples, const bool include, const int desired_samples) {
  for (const auto& record : SingleVariantReader{filename, samples, include})
    BOOST_CHECK_EQUAL(record.n_samples(), static_cast<uint32_t>(desired_samples));
}

BOOST_AUTO_TEST_CASE( single_variant_reader_sites_only )  
{
  single_variant_reader_sample_test("testdata/test_variants.vcf", vector<string>{}, true, 0); // exclude all samples (sites-only)
  // Fails: Issue #209
  //single_variant_reader_sample_test("testdata/test_variants.bcf", vector<string>{}, true, 0); // exclude all samples (sites-only)
}

BOOST_AUTO_TEST_CASE( single_variant_reader_include_all_samples ) 
{
  single_variant_reader_sample_test("testdata/test_variants.vcf", vector<string>{}, false, 3);  // include all samples by setting include == false and passing an empty list
  single_variant_reader_sample_test("testdata/test_variants.bcf", vector<string>{}, false, 3);  // include all samples by setting include == false and passing an empty list
}

BOOST_AUTO_TEST_CASE( single_variant_reader_including )  
{
  single_variant_reader_sample_test("testdata/test_variants.vcf", vector<string>{"NA12878"}, true, 1); // include only NA12878
  single_variant_reader_sample_test("testdata/test_variants.vcf", vector<string>{"NA12878", "NA12892"}, true, 2); // include both these samples
  single_variant_reader_sample_test("testdata/test_variants.bcf", vector<string>{"NA12878"}, true, 1); // include only NA12878
  single_variant_reader_sample_test("testdata/test_variants.bcf", vector<string>{"NA12878", "NA12892"}, true, 2); // include both these samples
}

BOOST_AUTO_TEST_CASE( single_variant_reader_excluding )  
{
  single_variant_reader_sample_test("testdata/test_variants.vcf", vector<string>{"NA12891"}, false, 2);  // exclude only NA12891
  single_variant_reader_sample_test("testdata/test_variants.vcf", vector<string>{"NA12891", "NA12878"}, false, 1);  // exclude both these samples
  single_variant_reader_sample_test("testdata/test_variants.bcf", vector<string>{"NA12891"}, false, 2);  // exclude only NA12891
  single_variant_reader_sample_test("testdata/test_variants.bcf", vector<string>{"NA12891", "NA12878"}, false, 1);  // exclude both these samples
}

BOOST_AUTO_TEST_CASE( single_variant_reader_missing_data )
{
  single_variant_reader_sample_test("testdata/test_variants_missing_data.vcf", vector<string>{"0001A","0002A","0003A","0004A","0005A","0006A","0007A"}, true, 7);       // include 7 samples
  single_variant_reader_sample_test("testdata/test_variants_missing_data.vcf", vector<string>{"0007B", "0008A", "0009A", "0009B"}, false, 7);   // exclude 4 samples

  for (const auto& record : SingleVariantReader{"testdata/test_variants_missing_data.vcf"}){ // include all 11 samples
    BOOST_CHECK_EQUAL(record.n_samples(), 11u); 
    const auto gt_for_all_samples = record.genotypes();
    const auto ad_for_all_samples = record.integer_individual_field("AD");
    const auto gq_for_all_samples = record.integer_individual_field("GQ");
    const auto pl_for_all_samples = record.integer_individual_field("PL");

    for (const auto& gt_for_single_sample: gt_for_all_samples) {
      BOOST_CHECK(missing(gt_for_single_sample));
      BOOST_CHECK(gt_for_single_sample.missing());
      BOOST_CHECK_EQUAL(gt_for_single_sample.size(), 2u);
      BOOST_CHECK(missing(gt_for_single_sample.allele_string(0)));
      BOOST_CHECK(missing(gt_for_single_sample.allele_string(1)));
      BOOST_CHECK(missing(gt_for_single_sample.allele_key(0)));
      BOOST_CHECK(missing(gt_for_single_sample.allele_key(1)));
      BOOST_CHECK(missing(gt_for_single_sample[0]));
      BOOST_CHECK(missing(gt_for_single_sample[1]));
    }
    for (const auto& ad_for_single_sample: ad_for_all_samples) {
      BOOST_CHECK(missing(ad_for_single_sample));
      BOOST_CHECK(ad_for_single_sample.missing());
    }
    for (const auto& gq_for_single_sample: gq_for_all_samples) {
      BOOST_CHECK(missing(gq_for_single_sample));
      BOOST_CHECK(gq_for_single_sample.missing());
    }
    for (const auto& pl_for_single_sample: pl_for_all_samples) {
      BOOST_CHECK(missing(pl_for_single_sample));
      BOOST_CHECK(pl_for_single_sample.missing());
    }
  }
}

BOOST_AUTO_TEST_CASE( single_variant_reader_vector )
{
  for (const auto& record : SingleVariantReader{vector<string>{"testdata/test_variants.vcf"}})
    BOOST_CHECK_EQUAL(record.n_samples(), 3u);
}

BOOST_AUTO_TEST_CASE( single_variant_reader_vector_all_samples )
{
  for (const auto& record : SingleVariantReader{vector<string>{"testdata/test_variants.vcf"}, vector<string>{}, false})  // include all samples by setting include == false and passing an empty list
    BOOST_CHECK_EQUAL(record.n_samples(), 3u);
}

BOOST_AUTO_TEST_CASE( single_variant_reader_vector_too_large )
{
  BOOST_CHECK_THROW((SingleVariantReader{vector<string>{"testdata/test_variants.vcf", "testdata/test_variants.vcf"}}), std::runtime_error);
}

BOOST_AUTO_TEST_CASE( variant_reader_move_constructor ) {
  auto r0 = SingleVariantReader{"testdata/test_variants.bcf"};
  auto r1 = SingleVariantReader{"testdata/test_variants.bcf"};
  auto m1 = check_move_constructor(r1);
  BOOST_CHECK_EQUAL(r0.begin().operator*().alignment_start(), m1.begin().operator*().alignment_start());
}

BOOST_AUTO_TEST_CASE( variant_iterator_move_test ) {
  auto reader0 = SingleVariantReader{"testdata/test_variants.vcf"};
  auto iter0 = reader0.begin();

  auto reader1 = SingleVariantReader{"testdata/test_variants.vcf"};
  auto iter1 = reader1.begin();
  auto moved = check_move_constructor(iter1);

  auto record0 = *iter0;
  auto moved_record = *moved;
  auto truth_index = 0u;
  for (auto record : {record0, moved_record})
    check_all_apis(record, truth_index);
}

const auto gvcf_truth_ref               = vector<string>{"T", "C", "GG"};
const auto gvcf_truth_chromosome        = vector<uint32_t>{0, 0, 1};
const auto gvcf_truth_alignment_starts  = vector<uint32_t>{10000000, 20000000, 10001000};
const auto gvcf_truth_alignment_stops   = vector<uint32_t>{10000000, 20000123, 10001001};
const auto gvcf_truth_n_alleles         = vector<uint32_t>{2, 2, 2};
const auto gvcf_truth_id                = vector<string>{"db2342", ".", "rs837472"};

BOOST_AUTO_TEST_CASE( gvcf_test_single ) {
  auto truth_index = 0u;
  const auto reader = SingleVariantReader{"testdata/test.g.vcf"};
  for (const auto& record : reader) {
    BOOST_CHECK_EQUAL(record.ref(), gvcf_truth_ref[truth_index]);
    BOOST_CHECK_EQUAL(record.chromosome(), gvcf_truth_chromosome[truth_index]);
    BOOST_CHECK_EQUAL(record.alignment_start(), gvcf_truth_alignment_starts[truth_index]);
    BOOST_CHECK_EQUAL(record.alignment_stop(), gvcf_truth_alignment_stops[truth_index]);
    BOOST_CHECK_EQUAL(record.n_alleles(), gvcf_truth_n_alleles[truth_index]);
    BOOST_CHECK_EQUAL(record.n_samples(), 3u);
    BOOST_CHECK_EQUAL(record.id(), gvcf_truth_id[truth_index]);
    ++truth_index;
  }
}

BOOST_AUTO_TEST_CASE( variant_reader_nonexistent_file ) {
  BOOST_CHECK_THROW(SingleVariantReader{"foo/bar/nonexistent.vcf"}, FileOpenException);
}

// MultipleVariantReader / MultipleVariantIterator
// see also multiple_variant_reader_test

auto truth_file_indices_3x5 = vector<unordered_set<uint32_t>> {{0,1,2},{0,1,2},{0,1,2},{0,1,2},{0,1,2},{0,1,2},{0,1,2}};

BOOST_AUTO_TEST_CASE( multiple_variant_reader_test ) {
  auto truth_index = 0u;
  const auto reader = MultipleVariantReader<MultipleVariantIterator>{{"testdata/test_variants.vcf",
                                                                      "testdata/mvar_rdr/test_variants.bcf",
                                                                      "testdata/mvar_rdr/test_variants.vcf.gz"}, false};
  for (const auto& vec : reader) {
    auto expected_file_indices = truth_file_indices_3x5[truth_index];
    for (const auto& pair : vec) {
      const auto& record = pair.first;
      check_all_apis(record, truth_index);

      auto find_result = expected_file_indices.find(pair.second);
      BOOST_CHECK(find_result != expected_file_indices.end());
      expected_file_indices.erase(find_result);
    }
    BOOST_CHECK(expected_file_indices.empty());   // check that we've seen and erased all expected
    ++truth_index;
  }
}

BOOST_AUTO_TEST_CASE( multiple_variant_reader_move_test ) {
  auto reader0 = MultipleVariantReader<MultipleVariantIterator>{{"testdata/test_variants.vcf", "testdata/mvar_rdr/test_variants.bcf"}, false};;
  auto reader1 = MultipleVariantReader<MultipleVariantIterator>{{"testdata/test_variants.vcf", "testdata/mvar_rdr/test_variants.bcf"}, false};;
  auto moved = check_move_constructor(reader1);

  auto record0 = reader0.begin().operator*();
  auto moved_record = moved.begin().operator*();
  auto truth_index = 0u;
  for (const auto& vec : {record0, moved_record}) {
    for (const auto& pair : vec) {
      const auto& record = pair.first;
      check_all_apis(record, truth_index);
    }
  }
}

BOOST_AUTO_TEST_CASE( multiple_variant_iterator_move_test ) {
  auto reader0 = MultipleVariantReader<MultipleVariantIterator>{{"testdata/test_variants.vcf", "testdata/mvar_rdr/test_variants.bcf"}, false};
  auto iter0 = reader0.begin();

  auto reader1 = MultipleVariantReader<MultipleVariantIterator>{{"testdata/test_variants.vcf", "testdata/mvar_rdr/test_variants.bcf"}, false};
  auto iter1 = reader1.begin();
  auto moved = check_move_constructor(iter1);

  auto record0 = *iter0;
  auto moved_record = *iter1;
  auto truth_index = 0u;
  for (auto vec : {record0, moved_record}) {
    for (const auto& pair : vec) {
      const auto& record = pair.first;
      check_all_apis(record, truth_index);
    }
  }
}

// IndexedVariantReader / IndexedVariantIterator
// see also indexed_variant_reader_test

/*
  NOTE: Updated test_variants.bcf and var_idx directory via-
    bcftools view testdata/test_variants.vcf -o testdata/var_idx/test_variants.bcf -O b
    bcftools view testdata/test_variants.vcf -o testdata/var_idx/test_variants_csi.vcf.gz -O z
    cp testdata/var_idx/test_variants_csi.vcf.gz testdata/var_idx/test_variants_tabix.vcf.gz
    bcftools index testdata/var_idx/test_variants.bcf
    bcftools index testdata/var_idx/test_variants_csi.vcf.gz
    bcftools index testdata/var_idx/test_variants_tabix.vcf.gz -t
    cp testdata/var_idx/test_variants.bcf testdata/test_variants.bcf
 */
const auto indexed_variant_vcf_inputs = vector<string>{"testdata/var_idx/test_variants_csi.vcf.gz", "testdata/var_idx/test_variants_tabix.vcf.gz"};
const auto indexed_variant_bcf_inputs = vector<string>{"testdata/var_idx/test_variants.bcf"};

const auto indexed_variant_chrom_full = vector<string> {"1", "20", "22"};
const auto indexed_variant_bp_full = vector<string> {"1:10000000-10000000", "20:10001000-10001000", "20:10002000-10002000", "20:10003000-10003000", "22:10004000-10004000"};
const auto indexed_variant_chrom_partial = vector<string> {"1"};
const auto indexed_variant_bp_partial = vector<string> {"20:10001000-10001000"};

// empty interval lists and full interval lists
BOOST_AUTO_TEST_CASE( indexed_variant_reader_full_test ) {
  for (const auto filename : indexed_variant_bcf_inputs) {
    for (const auto intervals : {vector<string>{}, IndexedVariantIterator::all_intervals, indexed_variant_chrom_full, indexed_variant_bp_full}) {
      auto truth_index = 0u;
      const auto reader = IndexedVariantReader<IndexedVariantIterator>{filename, intervals};
      for (const auto& record : reader) {
        check_all_apis(record, truth_index);
        ++truth_index;
      }
      BOOST_CHECK_EQUAL(truth_index, 5u);
    }
  }
}

// SyncedVariantReader / SyncedVariantIterator
// see also synced_variant_reader_test

const auto synced_variant_vcf_inputs = indexed_variant_vcf_inputs;
const auto synced_variant_bcf_inputs = indexed_variant_bcf_inputs;

const auto synced_variant_chrom_full = boost::algorithm::join(indexed_variant_chrom_full, ",");
const auto synced_variant_bp_full = boost::algorithm::join(indexed_variant_bp_full, ",");
const auto synced_variant_chrom_partial = boost::algorithm::join(indexed_variant_chrom_partial, ",");
const auto synced_variant_bp_partial = boost::algorithm::join(indexed_variant_bp_partial, ",");

// full interval lists
BOOST_AUTO_TEST_CASE( synced_variant_reader_full_test ) {
  for (const auto intervals : {string{""}, synced_variant_chrom_full, synced_variant_bp_full}) {
    for (const auto input_files : {synced_variant_vcf_inputs, synced_variant_bcf_inputs}) {
      auto truth_index = 0u;
      const auto reader = SyncedVariantReader<SyncedVariantIterator>{input_files, intervals};
      for (const auto& vec : reader) {
        for (const auto& record : vec)
          check_all_apis(record, truth_index);
        ++truth_index;
      }
      BOOST_CHECK_EQUAL(truth_index, 5u);
    }
  }
}
