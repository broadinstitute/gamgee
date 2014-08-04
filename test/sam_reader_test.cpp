#include "sam_reader.h"

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

BOOST_AUTO_TEST_CASE( move_constructor_and_assignment ) {
  auto reader1 = SingleSamReader{"testdata/test_simple.bam"};
  auto it1 = reader1.begin();
  auto reader2 = std::move(reader1);                                     // check move constructor
  reader1 = SingleSamReader{"testdata/test_simple.bam"};                 // check move assignment
  auto it2 = reader2.begin();
  auto it3 = reader1.begin();
  BOOST_CHECK_NE((*it1).alignment_start(), (*it2).alignment_start());    // these should be different because they are pointing at the exact same iterator, so it2 should be 1 record ahead it1
  BOOST_CHECK_EQUAL((*it1).alignment_start(), (*it3).alignment_start()); // these should be the same! both pointing at the first record
}


