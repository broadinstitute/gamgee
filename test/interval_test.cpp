#include <boost/test/unit_test.hpp>

#include "interval.h"
#include "test_utils.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

using namespace std;
using namespace gamgee;

void check_interval(const Interval& actual, const Interval& truth) {
  BOOST_CHECK_EQUAL(actual.chr(), truth.chr());
  BOOST_CHECK_EQUAL(actual.start(), truth.start());
  BOOST_CHECK_EQUAL(actual.stop(), truth.stop());
}

void check_interval_file(const string& filename, const uint32_t true_size, const Interval& true_interval)
{
  const auto intervals = vector<Interval>{read_intervals(filename)};
  BOOST_CHECK_EQUAL(intervals.size(), true_size);
  for (const auto& loc : intervals) {
    check_interval(loc, true_interval);
  }
}

BOOST_AUTO_TEST_CASE( read_interval_files ) 
{
  // truth data for all files (if files change, this has to change too)
  const auto true_size = 18;
  const auto true_interval = Interval{"20", 132342, 4832374};

  // test all files
  const auto test_files = vector<string>{"testdata/test_gatk.intervals", "testdata/test_bed.bed", "testdata/test_picard.interval_list", "testdata/test_mixed.loc"};
  for (auto& filename : test_files) 
    check_interval_file(filename, true_size, true_interval);

  // test singletons
  const auto singletons = vector<Interval>{read_intervals("testdata/test_singletons.intervals")};
  check_interval(singletons[0], Interval{"1",200,300});
  check_interval(singletons[1], Interval{"1",500,500});
  check_interval(singletons[2], Interval{"1",500,500});
}


void check_interval_conversion(const Interval& interval, const string& truth_file) {
  ostringstream test;
  ifstream truth_stream{truth_file}; 
  auto truth = string{}; 
  test << interval;
  while (truth_stream.peek() == '@')
    getline(truth_stream, truth);
  getline(truth_stream, truth);
  BOOST_CHECK_EQUAL(test.str(), truth);
}

BOOST_AUTO_TEST_CASE( interval_conversion ) 
{
  auto intervals = vector<Interval>{read_intervals("testdata/test_gatk.intervals")};
  auto& interval = intervals[0];
  check_interval_conversion(interval, "testdata/test_gatk.intervals");        // testing GATK intervals (as it is the default output
  interval.set_output_type(gamgee::Interval::IntervalType::PICARD);
  check_interval_conversion(interval, "testdata/test_picard.interval_list");  // testing PICARD intervals
  interval.set_output_type(gamgee::Interval::IntervalType::BED);
  check_interval_conversion(interval, "testdata/test_bed.bed");               // testing BED intervals
}

BOOST_AUTO_TEST_CASE( interval_tile )
{
  const auto i = Interval{"1", 200, 250};
  for (auto len = 10u; len <= 200; len += 10) {
    for (auto spacing = 1u; spacing + len < i.size(); ++spacing) {
      const auto left = i.tile_left(len,spacing);
      const auto right = i.tile_right(len,spacing);
      const auto expected_num_tiles = uint64_t(floor(i.size()) + spacing) / (len+spacing);  
      BOOST_CHECK_EQUAL(left.size(), expected_num_tiles);
      BOOST_CHECK_EQUAL(right.size(), expected_num_tiles);
      BOOST_CHECK_EQUAL(left[0].start(), i.start());
      BOOST_CHECK_EQUAL(right[0].stop(), i.stop());
      BOOST_CHECK_LE(left[left.size()-1].stop(), i.stop());
      BOOST_CHECK_GE(right[right.size()-1].start(), i.start());
      // check spacing between tiled targets
      for (auto i = 0u; i < left.size() - 1; ++i) {
        BOOST_CHECK_EQUAL(left[i+1].start() - left[i].stop() - 1, spacing); 
        BOOST_CHECK_EQUAL(right[i].start() - right[i+1].stop() - 1, spacing); 
      }
    }
  }
}

BOOST_AUTO_TEST_CASE( interval_tsca_tiling )
{
  const uint32_t flanking {30};
  Interval target {"1", 200, 2000};
  for (uint32_t insert_size = 50; insert_size < 500; insert_size += 50) {
    for (uint32_t spacing = 0; spacing < insert_size/2; spacing += 2) {  // avoid testing odd spacing/insert_size to avoid off-by-one rounding of the positions.
      vector<Interval> v {target.tsca_tiling(spacing, insert_size, flanking)};
      for (uint32_t i = 0; i < v.size() - 1; ++i) {
        BOOST_CHECK_EQUAL(v[i].stop() - v[i+1].start() + 1, (insert_size - spacing)/2);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE( interval_equality ) 
{
  Interval i {"TEST", 234, 432}, j {"TEST", 234, 432}; 
  BOOST_CHECK(i == j);
  i.set_start(23);
  BOOST_CHECK(!(i == j));
  i.set_start(234);
  j.set_stop(443);
  BOOST_CHECK(!(i == j));
  i.set_stop(234);
  j.set_stop(234);
  BOOST_CHECK(i == j);
  i.set_chr("TAST");
  BOOST_CHECK(!(i == j));
}

BOOST_AUTO_TEST_CASE( interval_copy_and_move_constructors ) {
  auto i0 = Interval {"A", 1000, 2000};
  auto copies = check_copy_constructor(i0);
  auto c2 = get<1>(copies);
  BOOST_CHECK(i0 ==  get<0>(copies));
  BOOST_CHECK(i0 == c2);
  BOOST_CHECK(i0 == get<2>(copies));
  c2.set_start(1500);
  BOOST_CHECK(i0 != c2);
  BOOST_CHECK(c2 != get<0>(copies));
  BOOST_CHECK(c2 != get<1>(copies));
  BOOST_CHECK(c2 != get<2>(copies)); 
  auto m1 = check_move_constructor(get<0>(copies));
  BOOST_CHECK(i0 == m1);
}

BOOST_AUTO_TEST_CASE( string_output )
{
  const auto testfile = "testdata/test_gatk.intervals";
  auto intervals = vector<Interval>{read_intervals(testfile)};
  const auto test_str = intervals[0].str();

  ifstream truth_stream{testfile};
  auto truth = string{};
  while (truth_stream.peek() == '@')
    getline(truth_stream, truth);
  getline(truth_stream, truth);
  BOOST_CHECK_EQUAL(test_str, truth);
}


