#include <boost/test/unit_test.hpp>

#include "sam.h"
#include "sam_reader.h"
#include "sam_builder.h"
#include "missing.h"

#include "test_utils.h"

#include <vector>

using namespace std;
using namespace gamgee;


void check_first_read(Sam& record) {
  const auto chr = 5u;
  const auto aln = 1000u;
  const auto expected_cigar = "76M";
  const auto expected_cigar_size = 1u;
  const auto expected_cigar_element_length = 76u;
  const auto expected_cigar_element_operator = CigarOperator::M;
  const auto expected_bases = "ACCCTAACCCTAACCCTAACCCTAACCATAACCCTAAGACTAACCCTAAACCTAACCCTCATAATCGAAATACAAC";
  const vector<uint8_t> expected_quals = {33, 33, 33, 33, 34, 31, 34, 30, 32, 32, 33, 34, 33, 33, 27, 21, 18, 29, 28, 33, 31, 29, 10, 33, 24, 12, 24, 10, 8, 17, 33, 23, 11, 10, 31, 18, 17, 22, 33, 20, 32, 29, 24, 15, 7, 7, 29, 12, 10, 6, 6, 18, 30, 7, 14, 6, 6, 6, 32, 8, 7, 6, 6, 16, 24, 7, 6, 22, 13, 11, 9, 9, 4, 8, 18, 25};
  const auto expected_mapq = 0;
  const auto expected_isize = -130;
  record.set_chromosome(chr);
  BOOST_CHECK_EQUAL(record.chromosome(), chr);
  record.set_alignment_start(aln);
  BOOST_CHECK_EQUAL(record.alignment_start(), aln);
  record.set_mate_chromosome(chr);
  BOOST_CHECK_EQUAL(record.mate_chromosome(), chr);
  record.set_mate_alignment_start(aln);
  BOOST_CHECK_EQUAL(record.mate_alignment_start(), aln);
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
  BOOST_CHECK_EQUAL(record.mapping_qual(), expected_mapq);
  record.set_mapping_qual(20);
  BOOST_CHECK_EQUAL(record.mapping_qual(), 20);
  BOOST_CHECK_EQUAL(record.insert_size(), expected_isize);
  record.set_insert_size(400);
  BOOST_CHECK_EQUAL(record.insert_size(), 400);
}

BOOST_AUTO_TEST_CASE( sam_body_simple_members_by_reference ) {
  for (auto& record : SingleSamReader {"testdata/test_simple.bam"}) {
    check_first_read(record);
    break;
  }
}

BOOST_AUTO_TEST_CASE( sam_body_simple_members_by_copy ) {
  for (auto record : SingleSamReader {"testdata/test_simple.bam"}) {
    check_first_read(record);
    break;
  }
}

BOOST_AUTO_TEST_CASE( sam_input_vector ) {
  for (auto& record : SingleSamReader {vector<string>{"testdata/test_simple.bam"}}) {  // should be replaced by mock sam_body
    check_first_read(record);
    break;
  }
}

BOOST_AUTO_TEST_CASE( sam_body_flags ) 
{
  auto record = *(SingleSamReader {"testdata/test_simple.bam"}.begin());
  record.set_paired();
  BOOST_CHECK(record.paired());
  record.set_not_paired();
  BOOST_CHECK(!record.paired());
  record.set_unmapped();
  BOOST_CHECK(record.unmapped());
  record.set_not_unmapped();
  BOOST_CHECK(!record.unmapped());
  record.set_mate_unmapped();
  BOOST_CHECK(record.mate_unmapped());
  record.set_not_mate_unmapped();
  BOOST_CHECK(!record.mate_unmapped());
  record.set_reverse();
  BOOST_CHECK(record.reverse());
  record.set_not_reverse();
  BOOST_CHECK(!record.reverse());
  record.set_mate_reverse();
  BOOST_CHECK(record.mate_reverse());
  record.set_not_mate_reverse();
  BOOST_CHECK(!record.mate_reverse());
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
}

BOOST_AUTO_TEST_CASE( sam_in_place_base_quals_modification ) {
  auto read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto read_quals = read.base_quals();
  auto expected_quals = vector<uint8_t>{4, 33, 50, 42, 34, 31};

  read_quals[0] = 4;
  read_quals[2] = 50;
  read_quals[3] = 42;

  for ( auto i = 0u; i < expected_quals.size(); ++i ) {
    BOOST_CHECK_EQUAL(read_quals[i], expected_quals[i]);
  }
}

BOOST_AUTO_TEST_CASE( sam_in_place_bases_modification ) {
  auto read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto read_bases = read.bases();
  auto expected_bases = vector<Base>{ Base::N, Base::C, Base::G, Base::C, Base::T, Base::C, Base::A, Base::N, Base::N, Base::T, Base::T, Base::A, Base::A };

  // These test modification of the lower 4 bits only, upper 4 bits only, and both upper and lower bits
  read_bases.set_base(0, Base::N);
  read_bases.set_base(2, Base::G);
  read_bases.set_base(5, Base::C);
  read_bases.set_base(7, Base::N);
  read_bases.set_base(8, Base::N);
  read_bases.set_base(9, Base::T);

  for ( auto i = 0u; i < expected_bases.size(); ++i ) {
    BOOST_CHECK_EQUAL(int(read_bases[i]), int(expected_bases[i]));
  }
}

BOOST_AUTO_TEST_CASE( sam_in_place_cigar_modification ) {
  auto read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto read_cigar = read.cigar();

  read_cigar[0] = Cigar::make_cigar_element(30, CigarOperator::I);

  BOOST_CHECK_EQUAL(read_cigar[0], Cigar::make_cigar_element(30, CigarOperator::I));
}

BOOST_AUTO_TEST_CASE( sam_cigar_templated_copy_and_move_constructors ) {
  auto it = SingleSamReader{"testdata/test_simple.bam"}.begin();
  const auto c0 = (*it).cigar();
  const auto copies = check_copy_constructor(c0);
  auto c1 = get<0>(copies);
  auto c2 = get<1>(copies);
  auto c3 = get<2>(copies);
  BOOST_CHECK(c0 == c1); 
  BOOST_CHECK(c0 == c2);
  BOOST_CHECK(c0 == c3);
  c1[0] = Cigar::make_cigar_element(30, CigarOperator::I);
  BOOST_CHECK(c0 != c1);
  auto m0 = (*it).cigar();
  auto m1 = check_move_constructor(m0);
  auto m2 = (*it).cigar();
  BOOST_CHECK(m1 == m2);
  m1[0] = Cigar::make_cigar_element(20, CigarOperator::I);
  BOOST_CHECK(m1 == m2);  // the underlying object is the same and it still exists 
  BOOST_CHECK(m0 == m1);  // the underlying object is the same and it still exists (hasn't been destroyed, so they must still match)
  BOOST_CHECK(m1 != c1);  // check that modifying the moved doesn't affect the copied
}

BOOST_AUTO_TEST_CASE( invalid_cigar_access ) {
  auto read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto cigar = read.cigar();
  BOOST_CHECK_THROW(cigar[cigar.size()] = 3, out_of_range); 
  BOOST_CHECK_THROW(cigar[cigar.size()], out_of_range); 
  auto x = 2; // just a random variable to assign to so we can test the const version of operator[]
  BOOST_CHECK_THROW(x = cigar[std::numeric_limits<uint32_t>::max()], out_of_range); 
  x++; // using it otherwise compiler will complain
}

BOOST_AUTO_TEST_CASE( comparing_different_cigars ) {
  const auto header = SingleSamReader{"testdata/test_simple.bam"}.header();
  auto builder = SamBuilder{header};
  builder.set_name("bla").set_bases("actg").set_base_quals({4,24,3,3});
  const auto read1 = builder.set_cigar("4M").build();
  const auto read2 = builder.set_cigar("2M1I1M").build();
  BOOST_CHECK(read1.cigar() != read2.cigar());
}

BOOST_AUTO_TEST_CASE( sam_base_quals_templated_copy_and_move_constructors ) {
  auto it = SingleSamReader{"testdata/test_simple.bam"}.begin();
  auto c0 = (*it).base_quals();
  auto copies = check_copy_constructor(c0);
  auto c1 = get<0>(copies);
  auto c2 = get<1>(copies);
  auto c3 = get<2>(copies);
  BOOST_CHECK(c0 == c1); 
  BOOST_CHECK(c0 == c2);
  BOOST_CHECK(c0 == c3);
  c1[0] = 99;
  BOOST_CHECK(c0 != c1);
  auto m0 = (*it).base_quals();
  auto m1 = check_move_constructor(m0);
  auto m2 = (*it).base_quals();
  BOOST_CHECK(m1 == m2);
  m1[0] = 90;
  BOOST_CHECK(m1 == m2);  // the underlying object is the same and it still exists 
  BOOST_CHECK(m0 == m1);  // the underlying object is the same and it still exists (hasn't been destroyed, so they must still match)
  BOOST_CHECK(m1 != c2);  // check that modifying the moved doesn't affect the copied
}

BOOST_AUTO_TEST_CASE( comparing_different_base_quals ) {
  const auto header = SingleSamReader{"testdata/test_simple.bam"}.header();
  auto builder = SamBuilder{header};
  builder.set_name("bla").set_bases("actg").set_cigar("4M");
  const auto read1 = builder.set_base_quals({4,4,3,2}).build();
  const auto read2 = builder.set_bases("act").set_cigar("3M").set_base_quals({4,4,3}).build(); // subtly different...
  BOOST_CHECK(read1.base_quals() != read2.base_quals());
}

BOOST_AUTO_TEST_CASE( sam_read_bases_templated_copy_and_move_constructors ) {
  auto it = SingleSamReader{"testdata/test_simple.bam"}.begin();
  auto c0 = (*it).bases();
  auto copies = check_copy_constructor(c0);
  auto c1 = get<0>(copies);
  auto c2 = get<1>(copies);
  auto c3 = get<2>(copies);
  BOOST_CHECK(c0 == c1); 
  BOOST_CHECK(c0 == c2);
  BOOST_CHECK(c0 == c3);
  c1.set_base(0, Base::C);
  BOOST_CHECK(c0 != c1);
  auto m0 = (*it).bases();
  auto m1 = check_move_constructor(m0);
  auto m2 = (*it).bases();
  BOOST_CHECK(m1 == m2);
  m1.set_base(0, Base::N);
  BOOST_CHECK(m1 == m2);  // the underlying object is the same and it still exists 
  BOOST_CHECK(m0 == m1);  // the underlying object is the same and it still exists (hasn't been destroyed, so they must still match)
  BOOST_CHECK(m1 != c2);  // check that modifying the moved doesn't affect the copied
}

BOOST_AUTO_TEST_CASE( comparing_different_read_bases ) {
  const auto header = SingleSamReader{"testdata/test_simple.bam"}.header();
  auto builder = SamBuilder{header};
  builder.set_name("bla").set_base_quals({4,4,3,2}).set_cigar("4M");
  const auto read1 = builder.set_bases("acgt").build();
  const auto read2 = builder.set_bases("act").set_cigar("3M").set_base_quals({4,4,3}).build(); // subtly different...
  BOOST_CHECK(read1.bases() != read2.bases());
}

BOOST_AUTO_TEST_CASE( invalid_read_bases_access ) {
  auto read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto bases = read.bases();
  BOOST_CHECK_THROW(bases[bases.size()], out_of_range); 
  BOOST_CHECK_THROW(bases[std::numeric_limits<uint32_t>::max()], out_of_range); 
}

BOOST_AUTO_TEST_CASE( sam_read_tags ) {
  const auto read1 = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  const auto read2 = *(SingleSamReader{"testdata/test_paired.bam"}.begin());

  const auto read1_pg_tag = read1.string_tag("PG");
  BOOST_CHECK_EQUAL(read1_pg_tag.name(), "PG");
  BOOST_CHECK_EQUAL(read1_pg_tag.value(), "0");
  BOOST_CHECK_EQUAL(read1_pg_tag.missing(), false);

  const auto read1_rg_tag = read1.string_tag("RG");
  BOOST_CHECK_EQUAL(read1_rg_tag.name(), "RG");
  BOOST_CHECK_EQUAL(read1_rg_tag.value(), "exampleBAM.bam");
  BOOST_CHECK_EQUAL(read1_rg_tag.missing(), false);

  const auto read1_sm_tag = read1.string_tag("SM");
  BOOST_CHECK_EQUAL(read1_sm_tag.name(), "SM");
  BOOST_CHECK_EQUAL(read1_sm_tag.value(), "exampleBAM.bam");
  BOOST_CHECK_EQUAL(read1_sm_tag.missing(), false);

  const auto read1_nonexistent_tag = read1.integer_tag("DR");
  BOOST_CHECK_EQUAL(read1_nonexistent_tag.missing(), true);

  const auto read2_nm_tag = read2.integer_tag("NM");
  BOOST_CHECK_EQUAL(read2_nm_tag.name(), "NM");
  BOOST_CHECK_EQUAL(read2_nm_tag.value(), 0);
  BOOST_CHECK_EQUAL(read2_nm_tag.missing(), false);

  const auto read2_md_tag = read2.string_tag("MD");
  BOOST_CHECK_EQUAL(read2_md_tag.name(), "MD");
  BOOST_CHECK_EQUAL(read2_md_tag.value(), "76");
  BOOST_CHECK_EQUAL(read2_md_tag.missing(), false);

  const auto read2_as_tag = read2.integer_tag("AS");
  BOOST_CHECK_EQUAL(read2_as_tag.name(), "AS");
  BOOST_CHECK_EQUAL(read2_as_tag.value(), 76);
  BOOST_CHECK_EQUAL(read2_as_tag.missing(), false);

  const auto read2_xs_tag = read2.integer_tag("XS");
  BOOST_CHECK_EQUAL(read2_xs_tag.name(), "XS");
  BOOST_CHECK_EQUAL(read2_xs_tag.value(), 0);
  BOOST_CHECK_EQUAL(read2_xs_tag.missing(), false);

  // check integer, char and float tags
  const auto read1_za_tag = read1.double_tag("ZA");
  BOOST_CHECK_EQUAL(read1_za_tag.name(), "ZA");
  BOOST_CHECK_CLOSE(read1_za_tag.value(), 2.3, 0.001);
  BOOST_CHECK_EQUAL(read1_za_tag.missing(), false);
  const auto read1_zb_tag = read1.integer_tag("ZB");
  BOOST_CHECK_EQUAL(read1_zb_tag.name(), "ZB");
  BOOST_CHECK_EQUAL(read1_zb_tag.value(), 23);
  BOOST_CHECK_EQUAL(read1_zb_tag.missing(), false);
  const auto read1_zc_tag = read1.char_tag("ZC");
  BOOST_CHECK_EQUAL(read1_zc_tag.name(), "ZC");
  BOOST_CHECK_EQUAL(read1_zc_tag.value(), 't');
  BOOST_CHECK_EQUAL(read1_zc_tag.missing(), false);

  // check missing functionality on missing tags
  const auto read2_nonexistent_string_tag = read2.string_tag("PP");
  BOOST_CHECK_EQUAL(read2_nonexistent_string_tag.missing(), true);
  BOOST_CHECK(missing(read2_nonexistent_string_tag));
  const auto read2_nonexistent_integer_tag = read2.integer_tag("PP");
  BOOST_CHECK_EQUAL(read2_nonexistent_integer_tag.missing(), true);
  BOOST_CHECK(missing(read2_nonexistent_integer_tag));
  const auto read2_nonexistent_double_tag = read2.double_tag("PP");
  BOOST_CHECK_EQUAL(read2_nonexistent_double_tag.missing(), true);
  BOOST_CHECK(missing(read2_nonexistent_double_tag));
  const auto read2_nonexistent_char_tag = read2.char_tag("PP");
  BOOST_CHECK_EQUAL(read2_nonexistent_char_tag.missing(), true);
  BOOST_CHECK(missing(read2_nonexistent_string_tag));

  // missing value due to type mismatches
  const auto not_a_char_tag = read1.char_tag("ZB");     // ZB is an integer tag
  BOOST_CHECK(missing(not_a_char_tag));              // this should yield "not a char" which is equal to a missing value
  const auto not_a_string_tag = read1.string_tag("ZB"); // ZB is an integer tag
  BOOST_CHECK(missing(not_a_string_tag));            // this should yield "not a char" which is equal to a missing value
}

BOOST_AUTO_TEST_CASE( sam_templated_copy_and_move_constructors ) {
  auto it = SingleSamReader{"testdata/test_simple.bam"}.begin();
  auto c0 = *it;
  auto copies = check_copy_constructor(c0);
  auto c1 = get<0>(copies);
  auto c2 = get<1>(copies);
  auto c3 = get<2>(copies);
  BOOST_CHECK_EQUAL(c0.alignment_start(), c1.alignment_start()); 
  BOOST_CHECK_EQUAL(c0.alignment_start(), c2.alignment_start());
  BOOST_CHECK_EQUAL(c0.alignment_start(), c3.alignment_start());
  c1.set_alignment_start(1);
  BOOST_CHECK_NE(c0.alignment_start(), c1.alignment_start());
  auto m0 = *it;
  auto m1 = check_move_constructor(m0);
  BOOST_CHECK(m0.empty());
  auto m2 = *it;
  BOOST_CHECK_EQUAL(m1.alignment_start(), m2.alignment_start());
  m1.set_alignment_start(5000);
  BOOST_CHECK_NE(m1.alignment_start(), m2.alignment_start());  // the underlying object is the same and it still exists 
  BOOST_CHECK_NE(m1.alignment_start(), c2.alignment_start());  // check that modifying the moved doesn't affect the copied
}

void check_read_alignment_starts_and_stops(const Sam& read, const uint32_t astart, const uint32_t astop, const uint32_t ustart, const uint32_t ustop) {
  BOOST_CHECK_EQUAL(read.alignment_start(), astart);
  BOOST_CHECK_EQUAL(read.alignment_stop(), astop);
  BOOST_CHECK_EQUAL(read.unclipped_start(), ustart);
  BOOST_CHECK_EQUAL(read.unclipped_stop(), ustop);
}

BOOST_AUTO_TEST_CASE( sam_unclipped_start_and_stop ) {
  const auto header = SingleSamReader{"testdata/test_simple.bam"}.header();
  auto builder = SamBuilder{header, false};
  const auto sam1 = builder.set_chromosome(0).set_alignment_start(100).set_cigar("20M").build();  // these are not complete reads, they only work for testing purposes!
  const auto sam2 = builder.set_chromosome(0).set_alignment_start(100).set_cigar("5S15M").build(); 
  const auto sam3 = builder.set_chromosome(0).set_alignment_start(100).set_cigar("15M5S").build(); 
  const auto sam4 = builder.set_chromosome(0).set_alignment_start(100).set_cigar("5S10M5S").build(); 
  check_read_alignment_starts_and_stops(sam1, 100, 119, 100, 119);
  check_read_alignment_starts_and_stops(sam2, 100, 114, 95, 114);
  check_read_alignment_starts_and_stops(sam3, 100, 114, 100, 119);
  check_read_alignment_starts_and_stops(sam4, 100, 109, 95, 114);
}

BOOST_AUTO_TEST_CASE( sam_input_vector_too_large ) {
  BOOST_CHECK_THROW((SingleSamReader {vector<string>{"testdata/test_simple.bam","testdata/test_simple.bam"}}), std::runtime_error);
}

BOOST_AUTO_TEST_CASE( sam_mate_operations ) {
  const auto truth_mate_alignment_stop = vector<uint32_t>{264,207,264,10486,10448,22924};
  const auto truth_mate_unclipped_start = vector<uint32_t>{255,198,255,10477,10397,22924};
  const auto truth_mate_unclipped_stop = vector<uint32_t>{264,207,264,10486,10458,22924};
  auto i = 0u;
  for (const auto& record : SingleSamReader{"testdata/test_simple.sam"}) {
    const auto tag = record.string_tag("MC");
    if (!missing(tag)) {
      BOOST_CHECK_EQUAL(record.mate_alignment_stop(tag), truth_mate_alignment_stop[i]);
      BOOST_CHECK_EQUAL(record.mate_alignment_stop(), truth_mate_alignment_stop[i]);  // check the other overload as well
      BOOST_CHECK_EQUAL(record.mate_unclipped_start(tag), truth_mate_unclipped_start[i]);
      BOOST_CHECK_EQUAL(record.mate_unclipped_start(), truth_mate_unclipped_start[i]);  // check the other overload as well
      BOOST_CHECK_EQUAL(record.mate_unclipped_stop(tag), truth_mate_unclipped_stop[i]);
      BOOST_CHECK_EQUAL(record.mate_unclipped_stop(), truth_mate_unclipped_stop[i++]);  // check the other overload as well
    }
    else {
      BOOST_CHECK_THROW(record.mate_alignment_stop(), invalid_argument);
      BOOST_CHECK_EQUAL(record.mate_alignment_stop(tag), record.mate_alignment_start());
      BOOST_CHECK_THROW(record.mate_unclipped_start(), invalid_argument);
      BOOST_CHECK_EQUAL(record.mate_unclipped_start(tag), record.mate_alignment_start());
      BOOST_CHECK_THROW(record.mate_unclipped_stop(), invalid_argument);
      BOOST_CHECK_EQUAL(record.mate_unclipped_stop(tag), record.mate_alignment_start());
    }
  }
}

BOOST_AUTO_TEST_CASE( sam_off_by_one_uber_test ) {
  const auto header = SingleSamReader{"testdata/test_simple.bam"}.header();
  auto builder = SamBuilder{header};
  const auto record = builder.set_name("TEST").set_bases("A").set_cigar("1M").set_base_quals({20}).set_alignment_start(1).build();
  BOOST_CHECK_EQUAL(record.alignment_start(), 1u);
  BOOST_CHECK_EQUAL(record.alignment_stop(), 1u);
}
