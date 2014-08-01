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

