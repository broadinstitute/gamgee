#include "multiple_variant_reader.h"
#include "reference_block_splitting_variant_iterator.h"

#include "test_utils.h"

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
//TODO const auto truth_refs = vector<string>{"T", "T", "T", "T", "T", "61", "81", "T", "T", "T", "151", "201", "T", "401", "T", "T", "701", "GG"};

BOOST_AUTO_TEST_CASE( split_reference_blocks )
{
  auto reader = GVCFReader{test_files, false};
  auto position_counter = 0u;
  for (const auto& vec : reader) {
    for (const auto& record : vec) {
      BOOST_CHECK_EQUAL(record.chromosome(), truth_contigs[position_counter]);
      BOOST_CHECK_EQUAL(record.alignment_start(), truth_block_starts[position_counter]);
      BOOST_CHECK_EQUAL(record.alignment_stop(), truth_block_stops[position_counter]);
      //TODO BOOST_CHECK_EQUAL(record.ref(), truth_refs[position_counter]);
    }
    ++position_counter;
  }
  BOOST_CHECK_EQUAL(position_counter, 18u);
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
  for (const auto vec : {record0, moved_record}) {
    for (const auto record : vec) {
      BOOST_CHECK_EQUAL(record.chromosome(), truth_contigs[position_counter]);
      BOOST_CHECK_EQUAL(record.alignment_start(), truth_block_starts[position_counter]);
      BOOST_CHECK_EQUAL(record.alignment_stop(), truth_block_stops[position_counter]);
      //TODO BOOST_CHECK_EQUAL(record.ref(), truth_refs[position_counter]);
    }
  }
}

// test cases based on problematic regions from real 1000 Genomes data
// mainly involving long overlapping reference alleles

const auto p1_file = "testdata/ref_block/problem1.vcf";
const auto p1_block_starts = vector<uint32_t>{11140505, 11140610, 11140611, 11140620, 11140621, 11140623, 11140629, 11140768};
const auto p1_block_stops = vector<uint32_t>{11140609, 11140610, 11140619, 11140628, 11140628, 11140628, 11140658, 11141021};

BOOST_AUTO_TEST_CASE( reference_block_problem_region_1 ) {
  auto counter = 0u;
  auto reader = GVCFReader{{p1_file}};
  for (const auto& vec : reader) {
    BOOST_CHECK(vec.size() == 1u);
    for (const auto& record : vec) {
      BOOST_CHECK(record.alignment_start() == p1_block_starts[counter]);
      BOOST_CHECK(record.alignment_stop() == p1_block_stops[counter]);
      counter++;
    }
  }
  BOOST_CHECK(counter == 8u);
}

const auto p2_files = vector<string>{"testdata/ref_block/problem2_file1.vcf", "testdata/ref_block/problem2_file2.vcf"};
const auto p2_block_starts = vector<uint32_t>{3319304, 3319319, 3319404, 3319407, 3319410,
  3319413, 3319433, 3319447, 3319558, 3319565, 3319577, 3319599, 3319600, 3319601, 3319602};
const auto p2_block_stops = vector<uint32_t>{3319318, 3319403, 3319406, 3319409, 3319412,
  3319432, 3319446, 3319557, 3319564, 3319576, 3319598, 3319599, 3319600, 3319601, 3319614};

BOOST_AUTO_TEST_CASE( reference_block_problem_region_2 ) {
  auto counter = 0u;
  auto reader = GVCFReader{p2_files};
  for (const auto& vec : reader) {
    for (const auto& record : vec) {
      BOOST_CHECK(record.alignment_start() == p2_block_starts[counter]);
      if (record.alt().size() == 1)
        BOOST_CHECK(record.alignment_stop() == p2_block_stops[counter]);
    }
    counter++;
  }

  BOOST_CHECK(counter == p2_block_starts.size());
}

const auto p3_files = vector<string>{"testdata/ref_block/problem3_file1.vcf", "testdata/ref_block/problem3_file2.vcf"};
const auto p3_block_starts = vector<uint32_t>{3319304, 3319319, 3319404, 3319407, 3319410,
  3319413, 3319433, 3319447, 3319558, 3319565, 3319577, 3319599, 3319600};
const auto p3_block_stops = vector<uint32_t>{3319318, 3319403, 3319406, 3319409, 3319412,
  3319432, 3319446, 3319557, 3319564, 3319576, 3319598, 3319599, 3319600};

BOOST_AUTO_TEST_CASE( reference_block_problem_region_3 ) {
  auto counter = 0u;
  auto reader = GVCFReader{p3_files};
  for (const auto& vec : reader) {
    for (const auto& record : vec) {
      BOOST_CHECK(record.alignment_start() == p3_block_starts[counter]);
      if (record.alt().size() == 1)
        BOOST_CHECK(record.alignment_stop() == p3_block_stops[counter]);
    }
    counter++;
  }

  BOOST_CHECK(counter == p3_block_starts.size());
}
