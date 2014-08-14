#include "indexed_sam_reader.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;
using namespace gamgee::utils;

BOOST_AUTO_TEST_CASE( indexed_single_readers )
{
  for (const auto& filename : {"testdata/test_simple.bam"}) {
    auto read_counter = 0u;
    for (const auto& sam : IndexedSingleSamReader{filename, contig_values::entire_file_list}) {
      BOOST_CHECK_EQUAL(sam.name().substr(0, 15), "30PPJAAXX090125");
      BOOST_CHECK_EQUAL(sam.chromosome(), 0);
      ++read_counter;
    }
    BOOST_CHECK_EQUAL(read_counter, 33u);
  }
}

BOOST_AUTO_TEST_CASE( vector_iterator_similarity_to_indexed_single_readers ) {
  // If this one fails, why is the parallel test below supposed to pass? begin is not expected to reset the start location in indexed readers??
  auto vector1 = std::vector<std::int32_t>{1, 2, 3};
  auto it1 = vector1.begin();
  auto vector2 = std::move(vector1); // check move constructor
  vector1 = std::vector<std::int32_t>{1, 2, 3}; // check move assignment
  auto it2 = vector2.begin();
  auto it3 = vector1.begin();
  BOOST_CHECK_NE((*it1), (*it2));    // these should be different because they are pointing at the exact same iterator, so it2 should be 1 record ahead it1
  BOOST_CHECK_EQUAL((*it1), (*it3)); // these should be the same! both pointing at the first record
}

BOOST_AUTO_TEST_CASE( indexed_single_readers_move_constructor_and_assignment ) {
  auto reader1 = IndexedSingleSamReader{"testdata/test_simple.bam", contig_values::entire_file_list};
  auto it1 = reader1.begin();
  auto reader2 = std::move(reader1); // check move constructor
  reader1 = IndexedSingleSamReader{"testdata/test_simple.bam", contig_values::entire_file_list}; // check move assignment
  auto it2 = reader2.begin();
  auto it3 = reader1.begin();
  BOOST_CHECK_NE((*it1).alignment_start(), (*it2).alignment_start());    // these should be different because they are pointing at the exact same iterator, so it2 should be 1 record ahead it1
  BOOST_CHECK_EQUAL((*it1).alignment_start(), (*it3).alignment_start()); // these should be the same! both pointing at the first record
}
