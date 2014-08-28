#include "sam_reader.h"

#include "test_utils.h"

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
    BOOST_CHECK_EQUAL(secondary_alignments, 7u);
    BOOST_CHECK_EQUAL(read_counter, 44u);
  }
}

BOOST_AUTO_TEST_CASE( single_sam_reader_move_test ) {
  auto reader0 = SingleSamReader{"testdata/test_simple.bam"};
  auto reader1 = SingleSamReader{"testdata/test_simple.bam"};
  auto moved = check_move_constructor(reader1);

  auto record0 = reader0.begin().operator*();
  auto moved_record = moved.begin().operator*();
  BOOST_CHECK_EQUAL(record0.name(), moved_record.name());
  BOOST_CHECK_EQUAL(record0.chromosome(), moved_record.chromosome());
}

BOOST_AUTO_TEST_CASE( sam_iterator_move_test ) {
  auto reader0 = SingleSamReader{"testdata/test_simple.bam"};
  auto iter0 = reader0.begin();

  auto reader1 = SingleSamReader{"testdata/test_simple.bam"};
  auto iter1 = reader1.begin();
  auto moved = check_move_constructor(iter1);

  auto record0 = *iter0;
  auto moved_record = *moved;
  BOOST_CHECK_EQUAL(record0.name(), moved_record.name());
  BOOST_CHECK_EQUAL(record0.chromosome(), moved_record.chromosome());
}
