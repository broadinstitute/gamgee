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

int count_records_vec(const string& filename) {
  auto counter = 0;
  for (auto& f : FastqReader{vector<string>{filename}}) {
    BOOST_CHECK_EQUAL(f.sequence(), "ACAAGAGATTTAAGAC");
    ++counter;
  }
  return counter;
}

void vector_too_large(const string& filename) {
  for (auto& f : FastqReader{vector<string>{filename, filename}}) {
    BOOST_CHECK_EQUAL(f.sequence(), "ACAAGAGATTTAAGAC");
  }
}

BOOST_AUTO_TEST_CASE( read_fastq ) 
{
  BOOST_CHECK_EQUAL(count_records("testdata/complete_same_seq.fa"), 2);
  BOOST_CHECK_EQUAL(count_records("testdata/complete_same_seq.fq"), 3);
}

BOOST_AUTO_TEST_CASE( read_fastq_vector )
{
  BOOST_CHECK_EQUAL(count_records_vec("testdata/complete_same_seq.fa"), 2);
  BOOST_CHECK_EQUAL(count_records_vec("testdata/complete_same_seq.fq"), 3);
}

BOOST_AUTO_TEST_CASE( read_fastq_vector_too_large )
{
  BOOST_CHECK_THROW(vector_too_large("testdata/complete_same_seq.fa"), std::runtime_error);
}
