#include "sam_reader.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( sam_header ) {
  auto reader = SingleSamReader{"testdata/test_simple.bam"};
  auto header = reader.header();
  BOOST_CHECK_EQUAL("chr1", header.sequence_name(0));
  BOOST_CHECK_EQUAL(100000, header.sequence_length("chr1"));
  BOOST_CHECK_EQUAL(100000, header.sequence_length(0));
  BOOST_CHECK_EQUAL(0, header.sequence_length("foo"));
  BOOST_CHECK_EQUAL(1, header.n_sequences());
}

// Need a way to modify the header in between these copies/moves to make sure these are working properly!
BOOST_AUTO_TEST_CASE( sam_header_constructors ) {
  auto reader = SingleSamReader{"testdata/test_simple.bam"};
  auto header1 = reader.header();
  auto header2 = header1;                  // copy constructor
  auto header3 = std::move(header1);       // move constructor
  header1 = header2;                       // copy assignment
  const auto header4 = std::move(header1); // transfer header1's memory somewhere else so we can reuse it for move assignment
  header1 = std::move(header3);            // move assignment
  header1 = header1;                       // self assignment
}
