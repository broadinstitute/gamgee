#include "../gamgee/fastq_reader.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

int count_records(const string& filename) {
  auto counter = 0;
  for (auto& f : FastqReader{filename}) {
    BOOST_CHECK_EQUAL(f.sequence(), "ACAAGAGATTTAAGAC");
    ++counter;
  }
  return counter;
}

BOOST_AUTO_TEST_CASE( read_fastq ) 
{
  BOOST_CHECK_EQUAL(count_records("testdata/complete_same_seq.fa"), 2);
  BOOST_CHECK_EQUAL(count_records("testdata/complete_same_seq.fq"), 3);
}
