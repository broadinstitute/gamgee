#include "reference_map.h"
#include "utils/utils.h"

#include <boost/test/unit_test.hpp>

#include <vector>
#include <string>
#include <unordered_map>

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( reference_map_constructor_test ) 
{
  // truth data for test reference (if file changes, this has to change too)
  const auto SEQUENCE = string{"AGGGTAGAGAGATAGAGATCCCCCCCCCCAGTACCNNNNAGTT"};
  const auto ref = ReferenceMap {"testdata/test_reference.fa"};
  for (auto& chr_seq : ref) {
    BOOST_CHECK_EQUAL(chr_seq.first.substr(0,3), "chr");
    BOOST_CHECK_EQUAL(chr_seq.second, SEQUENCE);
  }
}

BOOST_AUTO_TEST_CASE( get_sequence_test ) 
{
  // truth data for test reference (if file changes, this has to change too)
  const auto SEQUENCE = string{"AGGGTAGAGAGATAGAGATCCCCCCCCCCAGTACCNNNNAGTT"};
  const auto reference_map = ReferenceMap{"testdata/test_reference.fa"};
  for (auto start = 1u; start != SEQUENCE.length(); ++start) {
    for (auto len = 1u; len <= SEQUENCE.length() - start; ++len) {
      const auto interval = Interval{"chrA", start, start+len-1};
      BOOST_CHECK_EQUAL(reference_map.get_sequence(interval), SEQUENCE.substr(start-1, len));
      BOOST_CHECK_EQUAL(reference_map.get_sequence(interval, true), gamgee::utils::complement(SEQUENCE.substr(start-1, len)));
    }
  }
}

