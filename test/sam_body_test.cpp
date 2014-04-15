#include <boost/test/unit_test.hpp>
#include "../gamgee/sam_body.h"
#include "../gamgee/sam_reader.h"

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( sam_body_simple_members ) {
  const auto chr = 5u;
  const auto aln = 1000u;
  for (auto record : SingleSamReader {"testdata/test_simple.bam"}) {  // should be replaced by mock sam_body
    record.set_chromosome(chr);
    BOOST_CHECK_EQUAL(record.chromosome(), chr);
    record.set_alignment_start(aln);
    BOOST_CHECK_EQUAL(record.alignment_start(), aln);
    record.set_mate_chromosome(chr);
    BOOST_CHECK_EQUAL(record.mate_chromosome(), chr);
    record.set_mate_alignment_start(aln);
    BOOST_CHECK_EQUAL(record.mate_alignment_start(), aln);
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
