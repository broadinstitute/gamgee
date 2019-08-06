#include "variant/multiple_variant_reader.h"
#include "variant/reference_block_splitting_variant_iterator.h"

#include "test_utils.h"

#include <unordered_set>

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

using GVCFReader = MultipleVariantReader<ReferenceBlockSplittingVariantIterator>;

const auto test_files = vector<string>{"testdata/ref_block/test1.vcf", "testdata/ref_block/test2.vcf", "testdata/ref_block/test3.vcf",
  "testdata/ref_block/test4.vcf", "testdata/ref_block/test5.vcf"};
const auto truth_contigs = vector<uint32_t>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
const auto truth_block_starts = vector<uint32_t>{1, 20, 41, 42, 50, 61, 81, 101, 102, 103, 151, 201, 300, 401, 600, 650, 701, 10};
const auto truth_block_stops = vector<uint32_t>{19, 40, 41, 49, 60, 80, 100, 101, 102, 150, 200, 299, 400, 500, 649, 700, 750, 11};
// note: final block is not a ref block, but the block ends at start + 1 because the reference allele length is 2
const auto truth_refs = vector<string>{"A", "C", "G", "T", "A", "N", "N", "C", "G", "T", "N", "N", "A", "N", "C", "T", "N", "GG"};
auto truth_file_indices_split = vector<unordered_set<uint32_t>> {{0,2,3},{0,2,3,4},{0,2,3,4},{0,2,3,4},{0,1,2,3,4},{0,1,3,4},{0,1,3},{0,1,3},{0,1,3},{0,1,3},
        {0,3},{3},{1,3},{3},{4},{2,4},{2},{0,1,2,3,4}};

BOOST_AUTO_TEST_CASE( split_reference_blocks )
{
  auto reader = GVCFReader{test_files, false};
  auto truth_index = 0u;
  for (const auto& vec : reader) {
    auto expected_file_indices = truth_file_indices_split[truth_index];
    for (const auto& pair : vec) {
      const auto& record = pair.first;
      BOOST_CHECK_EQUAL(record.chromosome(), truth_contigs[truth_index]);
      BOOST_CHECK_EQUAL(record.alignment_start(), truth_block_starts[truth_index]);
      BOOST_CHECK_EQUAL(record.alignment_stop(), truth_block_stops[truth_index]);
      BOOST_CHECK_EQUAL(record.ref(), truth_refs[truth_index]);

      auto find_result = expected_file_indices.find(pair.second);
      BOOST_CHECK(find_result != expected_file_indices.end());
      expected_file_indices.erase(find_result);
    }
    BOOST_CHECK(expected_file_indices.empty());   // check that we've seen and erased all expected
    ++truth_index;
  }
  BOOST_CHECK_EQUAL(truth_index, truth_contigs.size());
}

template<bool v1, bool v2>
void test_lut_base()
{
  utils::MergedVCFLUTBase<v1, v2> lut {};
  lut.resize_luts_if_needed(10u, 14u);
  //Check size
  BOOST_CHECK_EQUAL(lut.m_num_input_vcfs, 10u);
  BOOST_CHECK_EQUAL(lut.m_num_merged_fields, 14u);
  //Everything should be invalid/missing
  for(auto i=0u;i<lut.m_num_input_vcfs;++i)
    for(auto j=0u;j<lut.m_num_merged_fields;++j)
    {
      BOOST_CHECK(gamgee::missing(lut.get_input_idx_for_merged(i,j)));
      BOOST_CHECK(gamgee::missing(lut.get_merged_idx_for_input(i,j)));
    }
  lut.add_input_merged_idx_pair(5u, 4, 7);
  //Everything should be invalid/missing, except for the pair added above
  for(auto i=0u;i<lut.m_num_input_vcfs;++i)
    for(auto j=0u;j<lut.m_num_merged_fields;++j)
    {
      if(i != 5u && j != 7u)
      {
	BOOST_CHECK(gamgee::missing(lut.get_input_idx_for_merged(i,j)));
      }
      if(i != 5u && j != 4u)
      {
	BOOST_CHECK(gamgee::missing(lut.get_merged_idx_for_input(i,j)));
      }
    }
  BOOST_CHECK_EQUAL(lut.get_input_idx_for_merged(5u,7), 4);
  BOOST_CHECK_EQUAL(lut.get_merged_idx_for_input(5u,4), 7);
  //resize
  lut.resize_luts_if_needed(12u, 17u);
  //Check size
  BOOST_CHECK_EQUAL(lut.m_num_input_vcfs, 12u);
  BOOST_CHECK_EQUAL(lut.m_num_merged_fields, 17u);
  //Everything should be invalid/missing, except for the pair added above
  for(auto i=0u;i<lut.m_num_input_vcfs;++i)
    for(auto j=0u;j<lut.m_num_merged_fields;++j)
    {
      if(i != 5u && j != 7u)
      {
	BOOST_CHECK(gamgee::missing(lut.get_input_idx_for_merged(i,j)));
      }
      if(i != 5u && j != 4u)
      {
	BOOST_CHECK(gamgee::missing(lut.get_merged_idx_for_input(i,j)));
      }
    }
  BOOST_CHECK_EQUAL(lut.get_input_idx_for_merged(5u,7), 4);
  BOOST_CHECK_EQUAL(lut.get_merged_idx_for_input(5u,4), 7);
  //reset 1 direction of the LUT
  lut.reset_merged_idx_for_input(5u, 4);
  BOOST_CHECK(gamgee::missing(lut.get_merged_idx_for_input(5u, 4)));
  BOOST_CHECK_EQUAL(lut.get_input_idx_for_merged(5u,7), 4);
  //reset other direction of the LUT
  lut.reset_input_idx_for_merged(5u, 7);
  BOOST_CHECK(gamgee::missing(lut.get_merged_idx_for_input(5u, 4)));
  BOOST_CHECK(gamgee::missing(lut.get_input_idx_for_merged(5u,7)));
}

BOOST_AUTO_TEST_CASE( merged_vcf_lut_test )
{
  test_lut_base<true, true>();
  test_lut_base<true, false>();
  test_lut_base<false, true>();
  test_lut_base<false, false>();
  utils::CombineAllelesLUT Z{5u};
  Z.add_input_merged_idx_pair(3u, 4, 8);
  BOOST_CHECK_EQUAL(Z.get_input_idx_for_merged(3u, 8), 4);
  BOOST_CHECK_EQUAL(Z.get_merged_idx_for_input(3u, 4), 8);
  BOOST_CHECK(gamgee::missing(Z.get_merged_idx_for_input(2u, 7)));
  Z.resize_luts_if_needed(15u);
  BOOST_CHECK_EQUAL(Z.get_input_idx_for_merged(3u, 8), 4);
  BOOST_CHECK_EQUAL(Z.get_merged_idx_for_input(3u, 4), 8);
  BOOST_CHECK(gamgee::missing(Z.get_merged_idx_for_input(2u, 7)));
}

enum VariantHeaderMergerTestPathEnum
{
  INPUT_FIELD_MISSING=0,
  INPUT_FIELD_VALID,
  INPUT_SAMPLE_MISSING,
  INPUT_SAMPLE_VALID,
  PL_ENUM_VALID,
  AD_ENUM_VALID,
  AF_ENUM_VALID,
  INVALID_ENUM_INVALID,
  TOTAL_NUM_VARIANT_HEADER_MERGER_TEST_PATHS
};

enum VariantHeaderMergerFieldsEnum
{
  PL_FIELD=0,
  AD_FIELD,
};

template<class VariantHeaderMergerTy>
void test_variant_header_merger(VariantHeaderMergerTy& merger, const VariantHeader& combined_header, const shared_ptr<bcf_hdr_t>& header,
    unsigned input_vcf_idx, vector<unsigned>& test_path_counters)
{
  //merged field-enum mappings
  merger.store_merged_field_idx_for_enum("PL", PL_FIELD);
  merger.store_merged_field_idx_for_enum("AD", AD_FIELD);
  //test large enum value
  merger.store_merged_field_idx_for_enum("AF", 100u);
  //Every field/sample in the merged header is either in the input header or isn't
  //Information in both the LUT and the header should match
  for(auto merged_field_idx=0u;merged_field_idx<1000u;++merged_field_idx)	//large number
  {
    for(auto field_type : {BCF_HL_FLT, BCF_HL_INFO, BCF_HL_FMT})
      if(combined_header.has_field(merged_field_idx, field_type))
      {
	auto input_field_idx = merger.get_input_header_idx_for_merged(input_vcf_idx, static_cast<int>(merged_field_idx));
	auto header_field_idx = bcf_hdr_id2int(header.get(), BCF_DT_ID, combined_header.get_field_name(merged_field_idx).c_str());
	if(gamgee::missing(input_field_idx))
	{
	  BOOST_CHECK_EQUAL(header_field_idx, -1);
	  ++test_path_counters[INPUT_FIELD_MISSING];
	}
	else
	{
	  BOOST_CHECK_EQUAL(header_field_idx, input_field_idx);
	  bool at_least_one_valid_type = false;
	  for(auto input_field_type : {BCF_HL_FLT, BCF_HL_INFO, BCF_HL_FMT})
	    at_least_one_valid_type = (at_least_one_valid_type ||  bcf_hdr_idinfo_exists(header.get(), input_field_type, input_field_idx));
	  BOOST_CHECK(at_least_one_valid_type);
	  BOOST_CHECK(bcf_hdr_int2id(header, BCF_DT_ID, input_field_idx));	//not NULL
	  BOOST_CHECK(bcf_hdr_int2id(header, BCF_DT_ID, input_field_idx) != "");
	  BOOST_CHECK_EQUAL(bcf_hdr_int2id(header, BCF_DT_ID, input_field_idx), combined_header.get_field_name(merged_field_idx));  
	  auto lut_merged_field_idx = merger.get_merged_header_idx_for_input(input_vcf_idx, input_field_idx);
	  BOOST_CHECK_EQUAL(lut_merged_field_idx, static_cast<int>(merged_field_idx));
	  ++test_path_counters[INPUT_FIELD_VALID];
	}
	auto enum_value = merger.get_enum_for_merged_field_idx(merged_field_idx);
	switch(enum_value)
	{
	  case PL_FIELD:
	    BOOST_CHECK_EQUAL(combined_header.get_field_name(merged_field_idx), "PL");
	    BOOST_CHECK_EQUAL(merged_field_idx, merger.get_merged_field_idx_for_enum(enum_value));
	    ++test_path_counters[PL_ENUM_VALID];
	    break;
	  case AD_FIELD:
	    BOOST_CHECK_EQUAL(combined_header.get_field_name(merged_field_idx), "AD");
	    BOOST_CHECK_EQUAL(merged_field_idx, merger.get_merged_field_idx_for_enum(enum_value));
	    ++test_path_counters[AD_ENUM_VALID];
	    break;
	  case 100:
	    BOOST_CHECK_EQUAL(combined_header.get_field_name(merged_field_idx), "AF");
	    BOOST_CHECK_EQUAL(merged_field_idx, merger.get_merged_field_idx_for_enum(enum_value));
	    ++test_path_counters[AF_ENUM_VALID];
	    break;
	  default:
	    const auto& field_name = combined_header.get_field_name(merged_field_idx);
	    BOOST_CHECK(field_name != "AF" && field_name != "AD" && field_name != "PL");
	    BOOST_CHECK_MESSAGE(enum_value == gamgee::missing_values::int32,
		"Valid enum mapping "<<enum_value<<" should not exist for field "<<field_name);
	    ++test_path_counters[INVALID_ENUM_INVALID];
	    break;
	}
      }
  }
  //Same check for samples
  for(auto merged_sample_idx=0u;merged_sample_idx<combined_header.n_samples();++merged_sample_idx)
  {
    auto input_sample_idx = merger.get_input_sample_idx_for_merged(input_vcf_idx, merged_sample_idx);
    auto header_sample_idx = bcf_hdr_id2int(header.get(), BCF_DT_SAMPLE, combined_header.get_sample_name(merged_sample_idx).c_str());
    if(gamgee::missing(input_sample_idx))
    {
      BOOST_CHECK_EQUAL(header_sample_idx, -1);
      ++test_path_counters[INPUT_SAMPLE_MISSING];
    }
    else
    {
      BOOST_CHECK_EQUAL(header_sample_idx, input_sample_idx);
      BOOST_CHECK(header->samples[input_sample_idx]);	//not null
      BOOST_CHECK(header->samples[input_sample_idx] != "");
      BOOST_CHECK_EQUAL(header->samples[input_sample_idx], combined_header.get_sample_name(merged_sample_idx));
      auto lut_merged_sample_idx = merger.get_merged_sample_idx_for_input(input_vcf_idx, input_sample_idx);
      BOOST_CHECK_EQUAL(lut_merged_sample_idx, merged_sample_idx);
      ++test_path_counters[INPUT_SAMPLE_VALID];
    }
  }
}

BOOST_AUTO_TEST_CASE( variant_header_merger_test )
{
  auto diff_hdr_test_files = vector<string>{ "testdata/ref_block/problem2_file2.vcf", "testdata/ref_block/test2.vcf" };
  auto many_hdr_test_files = vector<string>{
    "testdata/ref_block/problem2_file1.vcf",
    "testdata/mvr_hdr/test2.vcf",
    "testdata/var_hdr_merge/test1.vcf"
  };
  auto variant_header_merger_test_path_counters = vector<unsigned>(TOTAL_NUM_VARIANT_HEADER_MERGER_TEST_PATHS, 0u);
  for(auto const & multi_reader : { GVCFReader{test_files, false}, GVCFReader{diff_hdr_test_files, false}, GVCFReader{many_hdr_test_files, false} })
  {
    const vector<shared_ptr<bcf_hdr_t>>& header_vec = multi_reader.get_input_vcf_headers();
    const VariantHeader& combined_header = multi_reader.combined_header();
    //Duplicate hdr for testing FieldOrderedVariantHeaderMerger
    FieldOrderedVariantHeaderMerger field_ordered_merger;
    field_ordered_merger.add_headers(header_vec);
    VariantHeader duplicate_combined_header { field_ordered_merger.get_raw_merged_header() };

    InputOrderedVariantHeaderMerger& input_ordered_merger = 
      (const_cast<GVCFReader&>(multi_reader)).get_variant_header_merger();
    auto input_vcf_idx = 0u;
    for(auto& header : header_vec)
    {
      //Merged header should contain every sample and header field in input header
      for(auto i=0u;i<header->n[BCF_DT_ID];++i)
        for(auto field_type : {BCF_HL_FLT, BCF_HL_INFO, BCF_HL_FMT})
          if(bcf_hdr_idinfo_exists(header,field_type, i) && bcf_hdr_id2hrec(header,BCF_DT_ID,field_type,i))
              BOOST_CHECK(combined_header.has_field(bcf_hdr_int2id(header,BCF_DT_ID,i), field_type));
      for(auto i=0u;i<header->n[BCF_DT_SAMPLE];++i)
        if(header->samples[i] && bcf_hdr_id2int(header.get(), BCF_DT_SAMPLE, header->samples[i]) == i)
          BOOST_CHECK(combined_header.has_sample(header->samples[i]));
      test_variant_header_merger<InputOrderedVariantHeaderMerger>(
          input_ordered_merger, combined_header, header, input_vcf_idx,
          variant_header_merger_test_path_counters);
      test_variant_header_merger<FieldOrderedVariantHeaderMerger>(
          field_ordered_merger, duplicate_combined_header, header, input_vcf_idx,
          variant_header_merger_test_path_counters);
      ++input_vcf_idx;
    }
  }
  for(auto i=0u;i<variant_header_merger_test_path_counters.size();++i) {
    BOOST_CHECK_MESSAGE(variant_header_merger_test_path_counters[i] > 0u, "VariantHeaderMerger test path corresponding to "<<i<<" was not exercised\n");
  }
}

BOOST_AUTO_TEST_CASE( reference_block_iterator_move_test ) {
  auto reader0 = MultipleVariantReader<ReferenceBlockSplittingVariantIterator>{test_files, false};
  auto iter0 = reader0.begin();

  auto reader1 = MultipleVariantReader<ReferenceBlockSplittingVariantIterator>{test_files, false};
  auto iter1 = reader1.begin();
  auto moved = check_move_constructor(iter1);

  auto record0 = *iter0;
  auto moved_record = *iter1;
  const auto position_counter = 0;
  for (const auto& vec : {record0, moved_record}) {
    for (const auto& pair : vec) {
      const auto& record = pair.first;
      BOOST_CHECK_EQUAL(record.chromosome(), truth_contigs[position_counter]);
      BOOST_CHECK_EQUAL(record.alignment_start(), truth_block_starts[position_counter]);
      BOOST_CHECK_EQUAL(record.alignment_stop(), truth_block_stops[position_counter]);
      BOOST_CHECK_EQUAL(record.ref(), truth_refs[position_counter]);
    }
  }
}

// test cases based on problematic regions from real 1000 Genomes data
// mainly involving long overlapping reference alleles

const auto p1_file = "testdata/ref_block/problem1.vcf";
const auto p1_block_starts = vector<uint32_t>{11140505, 11140610, 11140611, 11140620, 11140621, 11140623, 11140629, 11140768};
const auto p1_block_stops = vector<uint32_t>{11140609, 11140610, 11140619, 11140628, 11140628, 11140628, 11140658, 11141021};

BOOST_AUTO_TEST_CASE( reference_block_problem_region_1 ) {
  auto truth_index = 0u;
  auto reader = GVCFReader{{p1_file}};
  for (const auto& vec : reader) {
    BOOST_CHECK(vec.size() == 1u);
    for (const auto& pair : vec) {
      const auto& record = pair.first;
      BOOST_CHECK(record.alignment_start() == p1_block_starts[truth_index]);
      BOOST_CHECK(record.alignment_stop() == p1_block_stops[truth_index]);
      truth_index++;
    }
  }
  BOOST_CHECK(truth_index == 8u);
}

const auto p2_files = vector<string>{"testdata/ref_block/problem2_file1.vcf", "testdata/ref_block/problem2_file2.vcf"};
const auto p2_block_starts = vector<uint32_t>{3319304, 3319319, 3319404, 3319407, 3319410,
  3319413, 3319433, 3319447, 3319558, 3319565, 3319577, 3319599, 3319600, 3319601, 3319602};
const auto p2_block_stops = vector<uint32_t>{3319318, 3319403, 3319406, 3319409, 3319412,
  3319432, 3319446, 3319557, 3319564, 3319576, 3319598, 3319599, 3319600, 3319601, 3319614};
auto truth_file_indices_p2 = vector<unordered_set<uint32_t>> {{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{1}};

BOOST_AUTO_TEST_CASE( reference_block_problem_region_2 ) {
  auto truth_index = 0u;
  auto reader = GVCFReader{p2_files};
  for (const auto& vec : reader) {
    auto expected_file_indices = truth_file_indices_p2[truth_index];
    for (const auto& pair : vec) {
      const auto& record = pair.first;
      BOOST_CHECK(record.alignment_start() == p2_block_starts[truth_index]);
      if (record.alt().size() == 1)
        BOOST_CHECK(record.alignment_stop() == p2_block_stops[truth_index]);
      auto find_result = expected_file_indices.find(pair.second);
      BOOST_CHECK(find_result != expected_file_indices.end());
      expected_file_indices.erase(find_result);
    }
    BOOST_CHECK(expected_file_indices.empty());   // check that we've seen and erased all expected
    truth_index++;
  }

  BOOST_CHECK(truth_index == p2_block_starts.size());
}

const auto p3_files = vector<string>{"testdata/ref_block/problem3_file1.vcf", "testdata/ref_block/problem3_file2.vcf"};
const auto p3_block_starts = vector<uint32_t>{3319304, 3319319, 3319404, 3319407, 3319410,
  3319413, 3319433, 3319447, 3319558, 3319565, 3319577, 3319599, 3319600};
const auto p3_block_stops = vector<uint32_t>{3319318, 3319403, 3319406, 3319409, 3319412,
  3319432, 3319446, 3319557, 3319564, 3319576, 3319598, 3319599, 3319600};
auto truth_file_indices_p3 = vector<unordered_set<uint32_t>> {{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0}};

BOOST_AUTO_TEST_CASE( reference_block_problem_region_3 ) {
  auto truth_index = 0u;
  auto reader = GVCFReader{p3_files};
  for (const auto& vec : reader) {
    auto expected_file_indices = truth_file_indices_p3[truth_index];
    for (const auto& pair : vec) {
      const auto& record = pair.first;
      BOOST_CHECK(record.alignment_start() == p3_block_starts[truth_index]);
      if (record.alt().size() == 1)
        BOOST_CHECK(record.alignment_stop() == p3_block_stops[truth_index]);

      auto find_result = expected_file_indices.find(pair.second);
      BOOST_CHECK(find_result != expected_file_indices.end());
      expected_file_indices.erase(find_result);
    }
    BOOST_CHECK(expected_file_indices.empty());   // check that we've seen and erased all expected
    truth_index++;
  }

  BOOST_CHECK(truth_index == p3_block_starts.size());
}

BOOST_AUTO_TEST_CASE( test_mvi_for_empty_vcf )
{
  auto mvr = MultipleVariantReader<MultipleVariantIterator> { vector<string>{"testdata/var_hdr_merge/empty_vcf.vcf"}, false };
  auto num_variants = 0u;
  for(auto it = mvr.begin(); it !=  mvr.end(); ++it)
    ++num_variants;
  BOOST_CHECK_EQUAL(num_variants, 0u);
  auto ref_split_mvr = GVCFReader { vector<string>{"testdata/var_hdr_merge/empty_vcf.vcf"}, false };
  for(auto it = ref_split_mvr.begin(); it !=  ref_split_mvr.end(); ++it)
    ++num_variants;
  BOOST_CHECK_EQUAL(num_variants, 0u);
}
