#include "sam_reader.h"
#include "sam_window_reader.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( single_readers ) 
{
  for (const auto& filename : {"testdata/test_simple.bam", "testdata/test_simple.sam"}) {
    auto read_counter = 0u;
    for (const auto& sam : SingleSamReader{filename}) {
      BOOST_CHECK_EQUAL(sam.name().substr(0, 15), "30PPJAAXX090125");
      BOOST_CHECK_EQUAL(sam.chromosome(), 0);
      ++read_counter;
    }
    BOOST_CHECK_EQUAL(read_counter, 33u);
  }
}

BOOST_AUTO_TEST_CASE( paired_readers ) {
  for (const auto& filename : {"testdata/test_paired.bam", "testdata/test_paired.sam"}) {
    auto read_counter = 0u;
    auto secondary_alignments = 0u;
    for (const auto& p : PairSamReader{filename}) {
      if (p.second.empty()) {
        BOOST_CHECK(p.first.secondary() || p.first.supplementary());
        ++secondary_alignments;
      }
      else {
        BOOST_CHECK(p.first.name() == p.second.name());
        read_counter += 2;
      }
    }
    BOOST_CHECK_EQUAL(secondary_alignments, 5u);
    BOOST_CHECK_EQUAL(read_counter, 42u);
  }
}

BOOST_AUTO_TEST_CASE( window_readers ) {
  const auto& filename = "testdata/test_simple.bam";
  auto record_count = 0;
  std::vector<int32_t> actual_record_counts = {};

  /*
  Full contract TBD.
  > Do windows overlap-- or end at-- start and end of contig?
  > Assuming we're 1-based, does this code have any off-by-one errors?
  */

  // TODO: Don't yet know syntax to effectively overload the same no-args begin() method signature with a different return value.
  // TODO: For now using the other copy/pasted SingleSamWindowReader.
  //for (const auto& sam_window : WindowSamReader{filename, 20000, 10000}) {
  for (const auto& sam_window : SingleSamWindowReader{filename, 20000, 10000}) {
    int32_t window_record_count = 0;
    for (const auto& sam : sam_window) {
      BOOST_CHECK_EQUAL(sam.chromosome(), 0);
      window_record_count++;
    }
    actual_record_counts.push_back(window_record_count);
    record_count += window_record_count;
  }

  std::vector<int32_t> expected_record_counts = {5, 7, 8, 7, 10, 8, 2, 4, 8, 4};
  BOOST_CHECK_EQUAL_COLLECTIONS(
      actual_record_counts.begin(), actual_record_counts.end(),
      expected_record_counts.begin(), expected_record_counts.end());
  BOOST_CHECK_EQUAL(record_count, 63);
}
