#include <boost/test/unit_test.hpp>

#include "fastq_reader.h"
#include "test_utils.h"

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

BOOST_AUTO_TEST_CASE( fastq_copy_and_move_constructor ) {
  auto it = FastqReader{"testdata/complete_same_seq.fa"}.begin();
  auto c0 = *it;
  auto copies = check_copy_constructor(c0);
  auto c1 = get<0>(copies);
  auto c2 = get<1>(copies);
  auto c3 = get<2>(copies);
  BOOST_CHECK(c0 == c1);
  BOOST_CHECK(c0 == c2);
  BOOST_CHECK(c0 == c3);
  c1.set_name("modified");
  BOOST_CHECK(c1 != c0);
  BOOST_CHECK(c1 != c2);
  auto m0 = *it;
  auto m1 = check_move_constructor(m0);
  auto m2 = *it;
  BOOST_CHECK(m1 == m2);
}
