#include <boost/test/unit_test.hpp>

#include "indexed_sam_reader.h"
#include "test_utils.h"


using namespace std;
using namespace gamgee;
using namespace gamgee::utils;

BOOST_AUTO_TEST_CASE( indexed_single_readers_intervals )
{
  for (const auto& filename : {"testdata/test_simple.bam"}) {
    auto read_counter = 0u;
    const auto interval_list = vector<string>{"chr1:201-257", "chr1:30001-40000", "chr1:59601-70000", "chr1:94001"};

    for (const auto& sam : IndexedSingleSamReader{filename, interval_list}) {
      BOOST_CHECK_EQUAL(sam.name().substr(0, 15), "30PPJAAXX090125");
      BOOST_CHECK_EQUAL(sam.chromosome(), 0);
      ++read_counter;
    }
    BOOST_CHECK_EQUAL(read_counter, 15u);

    const auto interval_read_counts = vector<int>{3, 4, 4, 4};
    auto interval_number = 0u;
    for (const auto& interval : interval_list) {
      auto interval_read_counter = 0u;
      for (const auto& sam : IndexedSingleSamReader{filename, vector<string>{interval}}) {
        BOOST_CHECK_EQUAL(sam.name().substr(0, 15), "30PPJAAXX090125");
        BOOST_CHECK_EQUAL(sam.chromosome(), 0);
        ++interval_read_counter;
      }
      BOOST_CHECK_EQUAL(interval_read_counter, interval_read_counts[interval_number]);
      ++interval_number;
    }
  }
}

BOOST_AUTO_TEST_CASE( indexed_single_readers_entire_file )
{
  const auto entire_file = vector<string>{"."};
  for (const auto& filename : {"testdata/test_simple.bam"}) {
    auto read_counter = 0u;
    for (const auto& sam : IndexedSingleSamReader{filename, entire_file}) {
      BOOST_CHECK_EQUAL(sam.name().substr(0, 15), "30PPJAAXX090125");
      BOOST_CHECK_EQUAL(sam.chromosome(), 0);
      ++read_counter;
    }
    BOOST_CHECK_EQUAL(read_counter, 33u);
  }
}

BOOST_AUTO_TEST_CASE( indexed_single_readers_unmapped )
{
  const auto unmapped_reads = vector<string>{"*"};
  for (const auto& filename : {"testdata/test_simple.bam"}) {
    auto read_counter = 0u;
    for (const auto& sam : IndexedSingleSamReader{filename, unmapped_reads}) {
      BOOST_CHECK_EQUAL(sam.name().substr(0, 15), "30PPJAAXX090125");
      BOOST_CHECK_EQUAL(sam.chromosome(), 0);
      ++read_counter;
    }
    BOOST_CHECK_EQUAL(read_counter, 0u);
  }
}

BOOST_AUTO_TEST_CASE( indexed_single_readers_empty )
{
  const auto empty_list = vector<string>{};
  for (const auto& filename : {"testdata/test_simple.bam"}) {
    auto read_counter = 0u;
    for (const auto& sam : IndexedSingleSamReader{filename, empty_list}) {
      BOOST_CHECK_EQUAL(sam.name().substr(0, 15), "30PPJAAXX090125");
      BOOST_CHECK_EQUAL(sam.chromosome(), 0);
      ++read_counter;
    }
    BOOST_CHECK_EQUAL(read_counter, 0u);
  }
}

BOOST_AUTO_TEST_CASE( indexed_single_readers_move_constructor_and_assignment ) {
  auto r0 = IndexedSingleSamReader{"testdata/test_simple.bam", vector<string>{"."}};
  auto m1 = check_move_constructor(r0);
  auto m2 = IndexedSingleSamReader{"testdata/test_simple.bam", vector<string>{"."}};
  BOOST_CHECK_EQUAL((*(m1.begin())).alignment_start(), (*(m2.begin())).alignment_start());
}

BOOST_AUTO_TEST_CASE( indexed_single_readers_begin_always_restarts ) {
  auto reader1 = IndexedSingleSamReader{"testdata/test_simple.bam", vector<string>{"."}};
  auto it1 = reader1.begin();
  ++it1;
  auto it2 = reader1.begin();
  BOOST_CHECK_NE((*it1).alignment_start(), (*it2).alignment_start());
  ++it2;
  BOOST_CHECK_EQUAL((*it1).alignment_start(), (*it2).alignment_start());
}

BOOST_AUTO_TEST_CASE( indexed_sam_iterator_move_test ) {
  auto reader0 = IndexedSingleSamReader{"testdata/test_simple.bam", vector<string>{"."}};
  auto iter0 = reader0.begin();

  auto reader1 = IndexedSingleSamReader{"testdata/test_simple.bam", vector<string>{"."}};
  auto iter1 = reader1.begin();
  auto moved = check_move_constructor(iter1);

  auto record0 = *iter0;
  auto moved_record = *moved;
  BOOST_CHECK_EQUAL(record0.name(), moved_record.name());
  BOOST_CHECK_EQUAL(record0.chromosome(), moved_record.chromosome());
}
