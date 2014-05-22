#include <boost/test/unit_test.hpp>
#include <vector>
#include "sam.h"
#include "sam_reader.h"

using namespace std;
using namespace gamgee;


BOOST_AUTO_TEST_CASE( sam_body_simple_members_by_reference ) {
  const auto chr = 5u;
  const auto aln = 1000u;
  const auto expected_cigar = "76M";
  const auto expected_cigar_size = 1;
  const auto expected_cigar_element_length = 76;
  const auto expected_cigar_element_operator = CigarOperator::M;
  const auto expected_bases = "ACCCTAACCCTAACCCTAACCCTAACCATAACCCTAAGACTAACCCTAAACCTAACCCTCATAATCGAAATACAAC";
  const vector<uint8_t> expected_quals = {33, 33, 33, 33, 34, 31, 34, 30, 32, 32, 33, 34, 33, 33, 27, 21, 18, 29, 28, 33, 31, 29, 10, 33, 24, 12, 24, 10, 8, 17, 33, 23, 11, 10, 31, 18, 17, 22, 33, 20, 32, 29, 24, 15, 7, 7, 29, 12, 10, 6, 6, 18, 30, 7, 14, 6, 6, 6, 32, 8, 7, 6, 6, 16, 24, 7, 6, 22, 13, 11, 9, 9, 4, 8, 18, 25};

  for (auto& record : SingleSamReader {"testdata/test_simple.bam"}) {  // should be replaced by mock sam_body
    record.set_chromosome(chr);
    BOOST_CHECK_EQUAL(record.chromosome(), chr);
    record.set_alignment_start(aln);
    BOOST_CHECK_EQUAL(record.alignment_start(), aln);
    record.set_mate_chromosome(chr);
    BOOST_CHECK_EQUAL(record.mate_chromosome(), chr);
    record.set_mate_alignment_start(aln);
    BOOST_CHECK_EQUAL(record.mate_alignment_start(), aln);

    // TODO: more comprehensive tests for variable-length data fields once setters are in place
    const auto actual_cigar = record.cigar();
    BOOST_CHECK_EQUAL(actual_cigar.to_string(), expected_cigar);
    BOOST_CHECK_EQUAL(actual_cigar.size(), expected_cigar_size);
    BOOST_CHECK_EQUAL(static_cast<int>(Cigar::cigar_op(actual_cigar[0])), static_cast<int>(expected_cigar_element_operator));
    BOOST_CHECK_EQUAL(Cigar::cigar_oplen(actual_cigar[0]), expected_cigar_element_length);

    BOOST_CHECK_EQUAL(record.bases().to_string(), expected_bases);

    const auto actual_quals = record.base_quals();
    BOOST_CHECK_EQUAL(actual_quals.size(), expected_quals.size());
    for ( auto i = 0u; i < actual_quals.size(); ++i ) {
      BOOST_CHECK_EQUAL(actual_quals[i], expected_quals[i]);
    }

    break;
  }
}

BOOST_AUTO_TEST_CASE( sam_body_simple_members_by_copy ) {
  const auto chr = 5u;
  const auto aln = 1000u;
  const auto expected_cigar = "76M";
  const auto expected_cigar_size = 1;
  const auto expected_cigar_element_length = 76;
  const auto expected_cigar_element_operator = CigarOperator::M;
  const auto expected_bases = "ACCCTAACCCTAACCCTAACCCTAACCATAACCCTAAGACTAACCCTAAACCTAACCCTCATAATCGAAATACAAC";
  const vector<uint8_t> expected_quals = {33, 33, 33, 33, 34, 31, 34, 30, 32, 32, 33, 34, 33, 33, 27, 21, 18, 29, 28, 33, 31, 29, 10, 33, 24, 12, 24, 10, 8, 17, 33, 23, 11, 10, 31, 18, 17, 22, 33, 20, 32, 29, 24, 15, 7, 7, 29, 12, 10, 6, 6, 18, 30, 7, 14, 6, 6, 6, 32, 8, 7, 6, 6, 16, 24, 7, 6, 22, 13, 11, 9, 9, 4, 8, 18, 25};

  for (auto record : SingleSamReader {"testdata/test_simple.bam"}) {  // should be replaced by mock sam_body
    record.set_chromosome(chr);
    BOOST_CHECK_EQUAL(record.chromosome(), chr);
    record.set_alignment_start(aln);
    BOOST_CHECK_EQUAL(record.alignment_start(), aln);
    record.set_mate_chromosome(chr);
    BOOST_CHECK_EQUAL(record.mate_chromosome(), chr);
    record.set_mate_alignment_start(aln);
    BOOST_CHECK_EQUAL(record.mate_alignment_start(), aln);

    // TODO: more comprehensive tests for variable-length data fields once setters are in place
    const auto actual_cigar = record.cigar();
    BOOST_CHECK_EQUAL(actual_cigar.to_string(), expected_cigar);
    BOOST_CHECK_EQUAL(actual_cigar.size(), expected_cigar_size);
    BOOST_CHECK_EQUAL(static_cast<int>(Cigar::cigar_op(actual_cigar[0])), static_cast<int>(expected_cigar_element_operator));
    BOOST_CHECK_EQUAL(Cigar::cigar_oplen(actual_cigar[0]), expected_cigar_element_length);

    BOOST_CHECK_EQUAL(record.bases().to_string(), expected_bases);

    const auto actual_quals = record.base_quals();
    BOOST_CHECK_EQUAL(actual_quals.size(), expected_quals.size());
    for ( auto i = 0u; i < actual_quals.size(); ++i ) {
      BOOST_CHECK_EQUAL(actual_quals[i], expected_quals[i]);
    }

    break;
  }
}

BOOST_AUTO_TEST_CASE( sam_body_flags ) 
{
  for (auto record : SingleSamReader {"testdata/test_simple.bam"}) {  // should be replaced by mock sam_body
    record.set_paired();
    BOOST_CHECK(record.paired());
    record.set_not_paired();
    BOOST_CHECK(!record.paired());
    record.set_unmapped();
    BOOST_CHECK(record.unmapped());
    record.set_not_unmapped();
    BOOST_CHECK(!record.unmapped());
    record.set_next_unmapped();
    BOOST_CHECK(record.next_unmapped());
    record.set_not_next_unmapped();
    BOOST_CHECK(!record.next_unmapped());
    record.set_reverse();
    BOOST_CHECK(record.reverse());
    record.set_not_reverse();
    BOOST_CHECK(!record.reverse());
    record.set_next_reverse();
    BOOST_CHECK(record.next_reverse());
    record.set_not_next_reverse();
    BOOST_CHECK(!record.next_reverse());
    record.set_first();
    BOOST_CHECK(record.first());
    record.set_not_first();
    BOOST_CHECK(!record.first());
    record.set_last();
    BOOST_CHECK(record.last());
    record.set_not_last();
    BOOST_CHECK(!record.last());
    record.set_secondary();
    BOOST_CHECK(record.secondary());
    record.set_not_secondary();
    BOOST_CHECK(!record.secondary());
    record.set_fail();
    BOOST_CHECK(record.fail());
    record.set_not_fail();
    BOOST_CHECK(!record.fail());
    record.set_duplicate();
    BOOST_CHECK(record.duplicate());
    record.set_not_duplicate();
    BOOST_CHECK(!record.duplicate());
    record.set_supplementary();
    BOOST_CHECK(record.supplementary());
    record.set_not_supplementary();
    BOOST_CHECK(!record.supplementary());
    break;  // only do this to one read
  }
}
