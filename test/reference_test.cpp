#include "reference_map.h"
#include "reference_iterator.h"

#include "utils/utils.h"

#include <boost/test/unit_test.hpp>

#include <vector>
#include <string>
#include <unordered_map>

using namespace std;
using namespace gamgee;

// truth data for test reference (if file changes, this has to change too)

const auto FILE1 = string{"testdata/test_reference.fa"};
const auto CHROMOSOMES1 = vector<string>{
  "chrA", "chrB", "chrC", "chrD", "chrE", "chrF", "chrG", "chrH",
  "chrI", "chrJ", "chrK", "chrL", "chrM", "chrN", "chrO", "chrQ", "chrR", "chrS", "chrT", "chrU"
};
const auto SEQ1 = string{"AGGGTAGAGAGATAGAGATCCCCCCCCCCAGTACCNNNNAGTT"};

const auto FILE2 = string{"testdata/test_reference2.fa"};
const auto CHROMOSOMES2 = vector<string>{ "chr1", "chr2" };
const auto CHR1_SEQ = string{"AGGGATCCCCCCCCCCAGTACCNNNNAGTT"};
const auto CHR2_SEQ = string{"NNNNNNAGGGATCCCNCCCCCCCAGTACCNNNNAGTT"};

BOOST_AUTO_TEST_CASE( reference_map_constructor_test ) 
{
  auto ref1 = ReferenceMap {FILE1};
  for (const auto& chr_seq : ref1) {
    BOOST_CHECK_EQUAL(chr_seq.first.substr(0,3), "chr");
    BOOST_CHECK_EQUAL(chr_seq.second, SEQ1);
  }

  auto ref2 = ReferenceMap {FILE2};
  BOOST_CHECK_EQUAL(ref2["chr1"], CHR1_SEQ);
  BOOST_CHECK_EQUAL(ref2["chr2"], CHR2_SEQ);
  // verify that it can access in arbitrary order
  BOOST_CHECK_EQUAL(ref2["chr1"], CHR1_SEQ);
}

BOOST_AUTO_TEST_CASE( reference_map_get_sequence_test )
{
  auto reference_map1 = ReferenceMap{FILE1};
  for (auto start = 1u; start != SEQ1.length(); ++start) {
    for (auto len = 1u; len <= SEQ1.length() - start; ++len) {
      const auto interval = Interval{"chrA", start, start+len-1};
      BOOST_CHECK_EQUAL(reference_map1.get_sequence(interval), SEQ1.substr(start-1, len));
      BOOST_CHECK_EQUAL(reference_map1.get_sequence(interval, true), gamgee::utils::complement(SEQ1.substr(start-1, len)));
    }
  }

  auto reference_map2 = ReferenceMap{FILE2};
  for (auto start = 1u; start != CHR1_SEQ.length(); ++start) {
    for (auto len = 1u; len <= CHR1_SEQ.length() - start; ++len) {
      const auto interval = Interval{"chr1", start, start+len-1};
      BOOST_CHECK_EQUAL(reference_map2.get_sequence(interval), CHR1_SEQ.substr(start-1, len));
      BOOST_CHECK_EQUAL(reference_map2.get_sequence(interval, true), gamgee::utils::complement(CHR1_SEQ.substr(start-1, len)));
    }
  }
  for (auto start = 1u; start != CHR2_SEQ.length(); ++start) {
    for (auto len = 1u; len <= CHR2_SEQ.length() - start; ++len) {
      const auto interval = Interval{"chr2", start, start+len-1};
      BOOST_CHECK_EQUAL(reference_map2.get_sequence(interval), CHR2_SEQ.substr(start-1, len));
      BOOST_CHECK_EQUAL(reference_map2.get_sequence(interval, true), gamgee::utils::complement(CHR2_SEQ.substr(start-1, len)));
    }
  }
}

BOOST_AUTO_TEST_CASE( reference_iterator_test ) {
  auto reference1 = ReferenceIterator{FILE1};
  for (const auto chr : CHROMOSOMES1) {
    for (auto counter = 1u; counter <= SEQ1.size(); counter++) {
      const char truth = SEQ1[counter - 1];
      BOOST_CHECK(truth == reference1.ref_base(chr, counter));
    }
  }

  auto reference2 = ReferenceIterator{FILE2};
  for (auto counter = 1u; counter <= CHR1_SEQ.size(); counter++) {
    const char truth = CHR1_SEQ[counter - 1];
    BOOST_CHECK(truth == reference2.ref_base("chr1", counter));
  }
  for (auto counter = 1u; counter <= CHR2_SEQ.size(); counter++) {
    const char truth = CHR2_SEQ[counter - 1];
    BOOST_CHECK(truth == reference2.ref_base("chr2", counter));
  }
}

