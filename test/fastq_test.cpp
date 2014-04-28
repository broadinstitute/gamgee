#include "fastq.h"
#include "fastq_reader.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp> 

#include <string>
#include <utility>

using namespace std;
using namespace gamgee;
using boost::test_tools::output_test_stream;

void check_fastq_output_with_file(const string& input, const string& truth) {
  output_test_stream fq_file{truth};       
  for (const auto& fq : FastqReader{input})
    fq_file << fq;
  BOOST_CHECK(fq_file.match_pattern());
}

void check_fastq_output(const string& name, const string& comment, const string& seq, const string& qual = "") {
  output_test_stream fq_output;
  auto fq = Fastq{name, comment, seq, qual};
  fq_output << fq;
  const auto fasta_truth = ">" + name + " " + comment + "\n" + seq + "\n";
  const auto fastq_truth = "@" + name + " " + comment + "\n" + seq + "\n+\n" + qual + "\n" ; 
  const auto truth = qual.empty() ?  move(fasta_truth) : move(fastq_truth); 
  BOOST_CHECK(fq_output.is_equal(truth));
}

void check_fastq_fields(const Fastq& record, const string& name, const string& comment, const string& sequence, const string& quals) {
  BOOST_CHECK_EQUAL(record.name()    , name);
  BOOST_CHECK_EQUAL(record.comment() , comment);  
  BOOST_CHECK_EQUAL(record.sequence(), sequence);
  BOOST_CHECK_EQUAL(record.quals()   , quals);
}

BOOST_AUTO_TEST_CASE( fastq_chop_test ) 
{
  const auto name    = string{"test"};
  const auto comment = string{"comm"};
  const auto seq     = string{"ACGTACGTACGT"};
  const auto qual    = string{"@#@$%$#@#$@$$#@!"};
  auto record = Fastq{name, comment, seq, qual};
  for (auto i = 0u; i != seq.length(); ++i) {
    record.chop(i);
    check_fastq_fields(record, name, comment, seq.substr(i), qual.substr(i));
    record = Fastq{name, comment, seq, qual};
  }
}

BOOST_AUTO_TEST_CASE( fastq_reverse_complement_test ) 
{
  const auto name    = string{"test"};
  const auto comment = string{"comm"};
  const auto seq     = string{"TTGATCTCCGAT"};
  const auto qual    = string{"@#@$%$#@#$@$$#@!"};
  const auto rev     = string{"ATCGGAGATCAA"};
  auto record = Fastq{name, comment, seq, qual};        // check that reversing only reverses the sequence (correctly)
  record.reverse_complement();
  check_fastq_fields(record, name, comment, rev, qual);
  record.reverse_complement();                          // re-reverse should go back to original
  BOOST_CHECK_EQUAL(record.sequence(), seq);
}

BOOST_AUTO_TEST_CASE( fastq_output )
{
  check_fastq_output("name", "comment", "ACAGAC", "!!$%#"); 
  check_fastq_output("name", "comment", "ACAGAC"); 
  check_fastq_output_with_file("testdata/test_clean.fq", "testdata/test_clean.fq"); // read a clean fastq
  check_fastq_output_with_file("testdata/complete_same_seq.fq", "testdata/test_clean.fq"); // read a dirty fastq
  check_fastq_output_with_file("testdata/complete_same_seq.fa", "testdata/complete_same_seq.fa"); // read a clean fasta
}

