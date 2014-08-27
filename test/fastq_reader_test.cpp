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

BOOST_AUTO_TEST_CASE( fastq_reader_move_constructor ) {
  auto reader0 = FastqReader{"testdata/complete_same_seq.fa"};
  auto reader1 = FastqReader{"testdata/complete_same_seq.fa"};
  auto moved = check_move_constructor(reader1);

  auto record0 = reader0.begin().operator*();
  auto moved_record = moved.begin().operator*();
  BOOST_CHECK(record0 == moved_record);
}

BOOST_AUTO_TEST_CASE( fastq_iterator_move_constructor ) {
  auto reader0 = FastqReader{"testdata/complete_same_seq.fa"};
  auto iter0 = reader0.begin();

  auto reader1 = FastqReader{"testdata/complete_same_seq.fa"};
  auto iter1 = reader1.begin();
  auto moved = check_move_constructor(iter1);

  auto record0 = *iter0;
  auto moved_record = *moved;
  BOOST_CHECK(record0 == moved_record);
}
