#include "sam_reader.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( sam_header )
{
	auto reader = SingleSamReader{"testdata/test_simple.bam"};
	auto header = reader.header();
  BOOST_CHECK_EQUAL("chr1", header.sequence_name(0));
  BOOST_CHECK_EQUAL(100000, header.sequence_length("chr1"));
  BOOST_CHECK_EQUAL(100000, header.sequence_length(0));
  BOOST_CHECK_EQUAL(0, header.sequence_length("foo"));
  BOOST_CHECK_EQUAL(1, header.n_sequences());
}
