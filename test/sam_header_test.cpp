#include <boost/test/unit_test.hpp>

#include "sam_reader.h"
#include "test_utils.h"

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( sam_header ) {
  auto reader = SingleSamReader{"testdata/test_simple.bam"};
  auto header = reader.header();
  BOOST_CHECK_EQUAL("chr1", header.sequence_name(0));
  BOOST_CHECK_EQUAL(100000u, header.sequence_length("chr1"));
  BOOST_CHECK_EQUAL(100000u, header.sequence_length(0));
  BOOST_CHECK_EQUAL(0u, header.sequence_length("foo"));
  BOOST_CHECK_EQUAL(1u, header.n_sequences());
}

BOOST_AUTO_TEST_CASE( sam_header_read_groups ) {
  auto reader = SingleSamReader{"testdata/test_simple.bam"};
  auto header = reader.header();
  auto rgs = header.read_groups();
  BOOST_CHECK_EQUAL(1u, rgs.size());
  BOOST_CHECK_EQUAL(rgs[0].id, "exampleBAM.bam");
  BOOST_CHECK_EQUAL(rgs[0].platform, "illumina");
  BOOST_CHECK_EQUAL(rgs[0].library, "exampleBAM.bam");
}

/** @todo Need a way to modify the header in between these copies/moves to make sure these are working properly! */
BOOST_AUTO_TEST_CASE( sam_header_constructors ) {
  auto reader = SingleSamReader{"testdata/test_simple.bam"};
  auto h0 = reader.header();
  auto copies = check_copy_constructor(h0);
  auto moves = check_copy_constructor(h0);
  // need builder to be able to modify the header and check. At least this test will blow up if something is not functional.
}
