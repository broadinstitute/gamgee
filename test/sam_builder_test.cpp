#include "sam/sam_builder.h"
#include "sam/sam_reader.h"

#include <boost/test/unit_test.hpp>
#include <limits>
#include<float.h>
#include <stdexcept>

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( set_name ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto original_read_name = original_read.name();
  auto builder = SamBuilder{original_read};
  auto modified_read = builder.set_name("new_name").build();

  // Verify that set_name worked, and that the original read is unaffected
  BOOST_CHECK_EQUAL(modified_read.name(), "new_name");
  BOOST_CHECK_EQUAL(original_read.name(), original_read_name);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.cigar() == original_read.cigar());
  BOOST_CHECK(modified_read.bases() == original_read.bases());
  BOOST_CHECK(modified_read.base_quals() == original_read.base_quals());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), original_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), original_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_cigar_by_vector ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto original_read_cigar = original_read.cigar().to_string();
  auto builder = SamBuilder{original_read};
  auto new_cigar = vector<CigarElement>{ Cigar::make_cigar_element(5, CigarOperator::M),
                                         Cigar::make_cigar_element(15, CigarOperator::I),
                                         Cigar::make_cigar_element(56, CigarOperator::M) };
  auto modified_read = builder.set_cigar(new_cigar).build();

  // Verify that set_cigar worked, and that the original read is unaffected
  BOOST_CHECK_EQUAL(modified_read.cigar().to_string(), "5M15I56M");
  BOOST_CHECK_EQUAL(original_read.cigar().to_string(), original_read_cigar);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == original_read.name());
  BOOST_CHECK(modified_read.bases() == original_read.bases());
  BOOST_CHECK(modified_read.base_quals() == original_read.base_quals());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), original_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), original_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_cigar_by_initializer_list ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto original_read_cigar = original_read.cigar().to_string();
  auto builder = SamBuilder{original_read};
  auto modified_read = builder.set_cigar({ Cigar::make_cigar_element(5, CigarOperator::M),
                                           Cigar::make_cigar_element(15, CigarOperator::I),
                                           Cigar::make_cigar_element(56, CigarOperator::M) }).build();

  // Verify that set_cigar worked, and that the original read is unaffected
  BOOST_CHECK_EQUAL(modified_read.cigar().to_string(), "5M15I56M");
  BOOST_CHECK_EQUAL(original_read.cigar().to_string(), original_read_cigar);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == original_read.name());
  BOOST_CHECK(modified_read.bases() == original_read.bases());
  BOOST_CHECK(modified_read.base_quals() == original_read.base_quals());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), original_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), original_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_cigar_by_string ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto original_read_cigar = original_read.cigar().to_string();
  auto builder = SamBuilder{original_read};
  auto modified_read = builder.set_cigar("5M15I56M").build();

  // Verify that set_cigar worked, and that the original read is unaffected
  BOOST_CHECK_EQUAL(modified_read.cigar().to_string(), "5M15I56M");
  BOOST_CHECK_EQUAL(original_read.cigar().to_string(), original_read_cigar);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == original_read.name());
  BOOST_CHECK(modified_read.bases() == original_read.bases());
  BOOST_CHECK(modified_read.base_quals() == original_read.base_quals());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), original_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), original_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_cigar_from_existing ) {
  auto starting_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto starting_read_cigar = starting_read.cigar().to_string();
  auto builder = SamBuilder{starting_read};
  auto modified_cigar = string{""};

  for ( auto cigar_donor : SingleSamReader{"testdata/test_paired.bam"} ) {
    auto cigar = cigar_donor.cigar();
    if ( cigar.to_string() != "76M" ) {
      builder.set_cigar(cigar);
      modified_cigar = cigar.to_string();
      break;
    }
  }
  auto modified_read = builder.build();

  // Verify that set_cigar worked, and that the original read is unaffected
  BOOST_CHECK_EQUAL(modified_read.cigar().to_string(), modified_cigar);
  BOOST_CHECK_EQUAL(starting_read.cigar().to_string(), starting_read_cigar);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == starting_read.name());
  BOOST_CHECK(modified_read.bases() == starting_read.bases());
  BOOST_CHECK(modified_read.base_quals() == starting_read.base_quals());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), starting_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), starting_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_empty_cigar ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto original_read_cigar = original_read.cigar().to_string();
  auto builder = SamBuilder{original_read};

  // Verify that set_cigar worked, for empty cigar string (e.g. in case of
  // unmapped read).
  auto modified_read_1 = builder.set_unmapped().set_cigar("").build();
  BOOST_CHECK_EQUAL(modified_read_1.cigar().to_string(), "");
  // Verify that original read is unaffected.
  BOOST_CHECK_EQUAL(original_read.cigar().to_string(), original_read_cigar);

  // Verify that set_cigar worked, for empty cigar string (e.g. in case of
  // unmapped read).
  auto modified_read_2 = builder.set_unmapped().set_cigar({}).build();
  BOOST_CHECK_EQUAL(modified_read_2.cigar().to_string(), "");
  // Verify that original read is unaffected.
  BOOST_CHECK_EQUAL(original_read.cigar().to_string(), original_read_cigar);
}

BOOST_AUTO_TEST_CASE( set_bases_by_vector ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto original_read_bases = original_read.bases().to_string();
  auto builder = SamBuilder{original_read};
  auto new_bases = vector<Base>{ Base::A, Base::C, Base::G, Base::T, Base::T, Base::C, Base::C, Base::G, Base::A, Base::A, Base::A, Base::C, Base::G, Base::T, Base::T, Base::C, Base::C, Base::G, Base::A, Base::A, Base::A, Base::C, Base::G, Base::T, Base::T, Base::C, Base::C, Base::G, Base::A, Base::A, Base::A, Base::C, Base::G, Base::T, Base::T, Base::C, Base::C, Base::G, Base::A, Base::A, Base::A, Base::C, Base::G, Base::T, Base::T, Base::C, Base::C, Base::G, Base::A, Base::A, Base::A, Base::C, Base::G, Base::T, Base::T, Base::C, Base::C, Base::G, Base::A, Base::A, Base::A, Base::C, Base::G, Base::T, Base::T, Base::C, Base::C, Base::G, Base::A, Base::A, Base::A, Base::C, Base::G, Base::T, Base::T, Base::C };
  auto modified_read = builder.set_bases(new_bases).build();

  // Verify that set_bases worked, and that the original read is unaffected
  BOOST_CHECK_EQUAL(modified_read.bases().to_string(), "ACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTC");
  BOOST_CHECK_EQUAL(original_read.bases().to_string(), original_read_bases);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == original_read.name());
  BOOST_CHECK(modified_read.cigar() == original_read.cigar());
  BOOST_CHECK(modified_read.base_quals() == original_read.base_quals());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), original_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), original_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_bases_by_initializer_list ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto original_read_bases = original_read.bases().to_string();
  auto builder = SamBuilder{original_read};
  auto modified_read = builder.set_bases({ Base::A, Base::C, Base::G, Base::T, Base::T, Base::C, Base::C, Base::G, Base::A, Base::A, Base::A, Base::C, Base::G, Base::T, Base::T, Base::C, Base::C, Base::G, Base::A, Base::A, Base::A, Base::C, Base::G, Base::T, Base::T, Base::C, Base::C, Base::G, Base::A, Base::A, Base::A, Base::C, Base::G, Base::T, Base::T, Base::C, Base::C, Base::G, Base::A, Base::A, Base::A, Base::C, Base::G, Base::T, Base::T, Base::C, Base::C, Base::G, Base::A, Base::A, Base::A, Base::C, Base::G, Base::T, Base::T, Base::C, Base::C, Base::G, Base::A, Base::A, Base::A, Base::C, Base::G, Base::T, Base::T, Base::C, Base::C, Base::G, Base::A, Base::A, Base::A, Base::C, Base::G, Base::T, Base::T, Base::C }).build();

  // Verify that set_bases worked, and that the original read is unaffected
  BOOST_CHECK_EQUAL(modified_read.bases().to_string(), "ACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTC");
  BOOST_CHECK_EQUAL(original_read.bases().to_string(), original_read_bases);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == original_read.name());
  BOOST_CHECK(modified_read.cigar() == original_read.cigar());
  BOOST_CHECK(modified_read.base_quals() == original_read.base_quals());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), original_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), original_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_bases_by_string ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto original_read_bases = original_read.bases().to_string();
  auto builder = SamBuilder{original_read};
  auto modified_read = builder.set_bases("ACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTC").build();

  // Verify that set_bases worked, and that the original read is unaffected
  BOOST_CHECK_EQUAL(modified_read.bases().to_string(), "ACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTCCGAAACGTTC");
  BOOST_CHECK_EQUAL(original_read.bases().to_string(), original_read_bases);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == original_read.name());
  BOOST_CHECK(modified_read.cigar() == original_read.cigar());
  BOOST_CHECK(modified_read.base_quals() == original_read.base_quals());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), original_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), original_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_bases_from_existing ) {
  auto reader = SingleSamReader{"testdata/test_simple.bam"};
  auto iter = reader.begin();
  auto starting_read = *iter;
  auto starting_read_bases = starting_read.bases().to_string();
  auto donor_read = ++iter;
  auto builder = SamBuilder{starting_read};
  auto modified_read = builder.set_bases(donor_read.bases()).build();

  // Verify that set_bases worked, and that the original read is unaffected
  BOOST_CHECK(modified_read.bases() == donor_read.bases());
  BOOST_CHECK_EQUAL(starting_read.bases().to_string(), starting_read_bases);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == starting_read.name());
  BOOST_CHECK(modified_read.cigar() == starting_read.cigar());
  BOOST_CHECK(modified_read.base_quals() == starting_read.base_quals());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), starting_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), starting_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_base_quals_by_vector ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto original_read_base_quals = original_read.base_quals().to_string();
  auto builder = SamBuilder{original_read};
  auto new_base_quals = vector<uint8_t>{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6 };
  auto modified_read = builder.set_base_quals(new_base_quals).build();

  // Verify that set_base_quals worked, and that the original read is unaffected
  auto modified_read_base_quals = modified_read.base_quals();
  for ( auto i = 0u; i < modified_read_base_quals.size(); ++i ) {
    BOOST_CHECK_EQUAL(modified_read_base_quals[i], new_base_quals[i]);
  }
  BOOST_CHECK_EQUAL(original_read.base_quals().to_string(), original_read_base_quals);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == original_read.name());
  BOOST_CHECK(modified_read.cigar() == original_read.cigar());
  BOOST_CHECK(modified_read.bases() == original_read.bases());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), original_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), original_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_base_quals_by_initializer_list ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto original_read_base_quals = original_read.base_quals().to_string();
  auto builder = SamBuilder{original_read};
  auto modified_read = builder.set_base_quals({ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6 }).build();

  // Verify that set_base_quals worked, and that the original read is unaffected
  BOOST_CHECK_EQUAL(modified_read.base_quals().to_string(), "1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6");
  BOOST_CHECK_EQUAL(original_read.base_quals().to_string(), original_read_base_quals);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == original_read.name());
  BOOST_CHECK(modified_read.cigar() == original_read.cigar());
  BOOST_CHECK(modified_read.bases() == original_read.bases());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), original_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), original_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_base_quals_from_existing ) {
  auto reader = SingleSamReader{"testdata/test_simple.bam"};
  auto iter = reader.begin();
  auto starting_read = *iter;
  auto starting_read_base_quals = starting_read.base_quals().to_string();
  auto donor_read = ++iter;
  auto builder = SamBuilder{starting_read};
  auto modified_read = builder.set_base_quals(donor_read.base_quals()).build();

  // Verify that set_base_quals worked, and that the original read is unaffected
  BOOST_CHECK(modified_read.base_quals() == donor_read.base_quals());
  BOOST_CHECK_EQUAL(starting_read.base_quals().to_string(), starting_read_base_quals);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == starting_read.name());
  BOOST_CHECK(modified_read.cigar() == starting_read.cigar());
  BOOST_CHECK(modified_read.bases() == starting_read.bases());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), starting_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), starting_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_multiple_data_fields ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto original_read_cigar = original_read.cigar().to_string();
  auto original_read_bases = original_read.bases().to_string();
  auto original_read_base_quals = original_read.base_quals().to_string();
  auto builder = SamBuilder{original_read};
  auto modified_read = builder.set_cigar("4M1I1M").set_bases("ACGTTA").set_base_quals({1, 2, 3, 4, 5, 6}).build();

  // Verify that the cigar, bases, and base quals were changed, and that the original read is unaffected
  BOOST_CHECK_EQUAL(modified_read.cigar().to_string(), "4M1I1M");
  BOOST_CHECK_EQUAL(modified_read.bases().to_string(), "ACGTTA");
  BOOST_CHECK_EQUAL(modified_read.base_quals().to_string(), "1 2 3 4 5 6");
  BOOST_CHECK_EQUAL(original_read.cigar().to_string(), original_read_cigar);
  BOOST_CHECK_EQUAL(original_read.bases().to_string(), original_read_bases);
  BOOST_CHECK_EQUAL(original_read.base_quals().to_string(), original_read_base_quals);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == original_read.name());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), original_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), original_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( add_tag ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto original_read_tags = original_read.all_tag_types();
  auto builder = SamBuilder{original_read};
  // construct numeric arrays
  std::vector<int64_t> int_array_values = {USHRT_MAX, 1, 2, 3 , 4};
  const auto numeric_uint16_array_tag =
    SamNumericArrayTag(SamTagType::UINTEGER16ARRAY, int_array_values);
  int_array_values.resize(5);
  for ( int i = 0; i < 5; i++ ) int_array_values[i] = -i;
  const auto numeric_int8_array_tag_2 =
    SamNumericArrayTag(SamTagType::INTEGER8ARRAY, int_array_values);
  std::vector<float> float_array_values = {-1, -2, -3 , -4, -5};
  const auto numeric_float_array_tag =
      SamNumericArrayTag(SamTagType::FLOATARRAY, float_array_values);

  auto modified_read = builder.clear_tags().
      add_char_tag("aa", 'C').
      add_integer_tag("bb", INT_MAX).
      add_integer_tag("cc", INT_MIN).
      add_integer_tag("dd", UINT_MAX).
      add_float_tag("ee", FLT_MAX).
      add_float_tag("ff", FLT_MIN).
      add_double_tag("gg", DBL_MAX).
      add_double_tag("hh", DBL_MIN).
      add_string_tag("ii", "string tag.").
      add_string_tag("jj", "second string tag.").
      // An example of byte array [0x1a, 0xe3, 0x1]
      add_byte_array_tag("kk", "1AE301").
      add_numeric_array_tag("ll", numeric_uint16_array_tag).
      add_numeric_array_tag("mm", numeric_int8_array_tag_2).
      add_numeric_array_tag("nn", numeric_float_array_tag).
      build();

  // Verify that add_*_tag worked, and that the original read is unaffected
  BOOST_CHECK_EQUAL(modified_read.all_tag_types().size(), 14);
  BOOST_CHECK_EQUAL(modified_read.char_tag("aa").value(), 'C');
  BOOST_CHECK_EQUAL(modified_read.integer_tag("bb").value(), INT_MAX);
  BOOST_CHECK_EQUAL(modified_read.integer_tag("cc").value(), INT_MIN);
  BOOST_CHECK_EQUAL(modified_read.integer_tag("dd").value(), UINT_MAX);
  BOOST_CHECK_EQUAL(modified_read.double_tag("ee").value(), FLT_MAX);
  BOOST_CHECK_EQUAL(modified_read.double_tag("ff").value(), FLT_MIN);
  BOOST_CHECK_EQUAL(modified_read.double_tag("gg").value(), DBL_MAX);
  BOOST_CHECK_EQUAL(modified_read.double_tag("hh").value(), DBL_MIN);
  BOOST_CHECK_EQUAL(modified_read.string_tag("ii").value(), "string tag.");
  BOOST_CHECK_EQUAL(modified_read.string_tag("jj").value(), "second string tag.");
  BOOST_CHECK_EQUAL(modified_read.byte_array_tag("kk").value(), "1AE301");
  BOOST_CHECK(modified_read.numeric_array_tag("ll").value() == numeric_uint16_array_tag);
  BOOST_CHECK(modified_read.numeric_array_tag("mm").value() == numeric_int8_array_tag_2);
  BOOST_CHECK(modified_read.numeric_array_tag("nn").value() == numeric_float_array_tag);
  BOOST_CHECK(original_read.all_tag_types() == original_read_tags);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == original_read.name());
  BOOST_CHECK(modified_read.bases() == original_read.bases());
  BOOST_CHECK(modified_read.base_quals() == original_read.base_quals());
  BOOST_CHECK(modified_read.cigar() == original_read.cigar());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), original_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), original_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_tag_from_existing ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  for ( auto read : SingleSamReader{"testdata/test_paired.bam"} ) {
    if (original_read.all_tag_types().size() > 0) break;
    original_read = read;
  }
  const auto original_read_tags = original_read.all_tag_types();
  // Test requires an alignemnt with aux tags.
  BOOST_CHECK(original_read_tags.size() != 0);
  auto builder = SamBuilder{original_read};
  auto modified_read = builder.add_integer_tag("aa", 1).build();

  // Verify that set_sam_tags worked, and that the original read is unaffected
  BOOST_CHECK_EQUAL(modified_read.integer_tag("aa").value(), 1);
  BOOST_CHECK_EQUAL(modified_read.all_tag_types().size(), original_read_tags.size() + 1);
  for (const auto& tag : original_read_tags) {
    switch (tag.second) {
      case SamTagType::INTEGER:
        BOOST_CHECK_EQUAL(modified_read.integer_tag(tag.first).value(),
                          original_read.integer_tag(tag.first).value());
        break;
      case SamTagType::DOUBLE:
        BOOST_CHECK_EQUAL(modified_read.double_tag(tag.first).value(),
                          original_read.double_tag(tag.first).value());
        break;
      case SamTagType::CHAR:
        BOOST_CHECK_EQUAL(modified_read.char_tag(tag.first).value(),
                          original_read.char_tag(tag.first).value());
        break;
      case SamTagType::STRING:
        BOOST_CHECK(modified_read.string_tag(tag.first).value() ==
                    original_read.string_tag(tag.first).value());
        break;
    }

  }

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == original_read.name());
  BOOST_CHECK(modified_read.bases() == original_read.bases());
  BOOST_CHECK(modified_read.base_quals() == original_read.base_quals());
  BOOST_CHECK(modified_read.cigar() == original_read.cigar());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), original_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), original_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_multiple_data_fields_one_time_build ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto original_read_cigar = original_read.cigar().to_string();
  auto original_read_bases = original_read.bases().to_string();
  auto original_read_base_quals = original_read.base_quals().to_string();
  auto builder = SamBuilder{original_read};
  builder.set_cigar("4M1I1M").set_bases("ACGTTA").set_base_quals({1, 2, 3, 4, 5, 6});
  auto modified_read = builder.one_time_build();

  // Verify that the cigar, bases, and base quals were changed, and that the original read is unaffected
  BOOST_CHECK_EQUAL(modified_read.cigar().to_string(), "4M1I1M");
  BOOST_CHECK_EQUAL(modified_read.bases().to_string(), "ACGTTA");
  BOOST_CHECK_EQUAL(modified_read.base_quals().to_string(), "1 2 3 4 5 6");
  BOOST_CHECK_EQUAL(original_read.cigar().to_string(), original_read_cigar);
  BOOST_CHECK_EQUAL(original_read.bases().to_string(), original_read_bases);
  BOOST_CHECK_EQUAL(original_read.base_quals().to_string(), original_read_base_quals);

  // Verify that other fields were inherited from the original read
  BOOST_CHECK(modified_read.name() == original_read.name());
  BOOST_CHECK_EQUAL(modified_read.chromosome(), original_read.chromosome());
  BOOST_CHECK_EQUAL(modified_read.alignment_start(), original_read.alignment_start());
}

BOOST_AUTO_TEST_CASE( set_multiple_data_fields_from_scratch ) {
  auto header = SingleSamReader{"testdata/test_simple.bam"}.header();
  auto builder = SamBuilder{header};
  builder.set_name("foo").set_cigar("1M1I1M").set_bases("ACT").set_base_quals({1, 2, 3});
  builder.set_chromosome(1).set_alignment_start(142).set_mate_alignment_start(199);
  auto constructed_read = builder.build();

  BOOST_CHECK_EQUAL(constructed_read.name(), "foo");
  BOOST_CHECK_EQUAL(constructed_read.cigar().to_string(), "1M1I1M");
  BOOST_CHECK_EQUAL(constructed_read.bases().to_string(), "ACT");
  BOOST_CHECK_EQUAL(constructed_read.base_quals().to_string(), "1 2 3");
  BOOST_CHECK_EQUAL(constructed_read.chromosome(), 1u);
  BOOST_CHECK_EQUAL(constructed_read.alignment_start(), 142u);
  BOOST_CHECK_EQUAL(constructed_read.mate_alignment_start(), 199u);
}

BOOST_AUTO_TEST_CASE( reconstruct_complete_read ) {
  auto original_read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());

  // Start from scratch with only the header, and see if we can construct a full copy of original_read
  auto builder = SamBuilder{original_read.header()};

  builder.set_name(original_read.name()).set_cigar(original_read.cigar()).set_bases(original_read.bases()).set_base_quals(original_read.base_quals());
  builder.set_chromosome(original_read.chromosome()).set_alignment_start(original_read.alignment_start());
  builder.set_mate_chromosome(original_read.mate_chromosome()).set_mate_alignment_start(original_read.mate_alignment_start());
  original_read.paired() ? builder.set_paired() : builder.set_not_paired();
  original_read.unmapped() ? builder.set_unmapped() : builder.set_not_unmapped();
  original_read.mate_unmapped() ? builder.set_mate_unmapped() : builder.set_not_mate_unmapped();
  original_read.reverse() ? builder.set_reverse() : builder.set_not_reverse();
  original_read.mate_reverse() ? builder.set_mate_reverse() : builder.set_not_mate_reverse();
  original_read.first() ? builder.set_first() : builder.set_not_first();
  original_read.last() ? builder.set_last() : builder.set_not_last();
  original_read.secondary() ? builder.set_secondary() : builder.set_not_secondary();
  original_read.fail() ? builder.set_fail() : builder.set_not_fail();
  original_read.duplicate() ? builder.set_duplicate() : builder.set_not_duplicate();
  original_read.supplementary() ? builder.set_supplementary() : builder.set_not_supplementary();

  auto read = builder.build();

  BOOST_CHECK_EQUAL(read.name(), original_read.name());
  BOOST_CHECK(read.cigar() == original_read.cigar());
  BOOST_CHECK(read.bases() == original_read.bases());
  BOOST_CHECK(read.base_quals() == original_read.base_quals());
  BOOST_CHECK_EQUAL(read.chromosome(), original_read.chromosome());
  BOOST_CHECK_EQUAL(read.alignment_start(), original_read.alignment_start());
  BOOST_CHECK_EQUAL(read.mate_chromosome(), original_read.mate_chromosome());
  BOOST_CHECK_EQUAL(read.mate_alignment_start(), original_read.mate_alignment_start());
  BOOST_CHECK_EQUAL(read.paired(), original_read.paired());
  BOOST_CHECK_EQUAL(read.unmapped(), original_read.unmapped());
  BOOST_CHECK_EQUAL(read.mate_unmapped(), original_read.mate_unmapped());
  BOOST_CHECK_EQUAL(read.reverse(), original_read.reverse());
  BOOST_CHECK_EQUAL(read.mate_reverse(), original_read.mate_reverse());
  BOOST_CHECK_EQUAL(read.first(), original_read.first());
  BOOST_CHECK_EQUAL(read.last(), original_read.last());
  BOOST_CHECK_EQUAL(read.secondary(), original_read.secondary());
  BOOST_CHECK_EQUAL(read.fail(), original_read.fail());
  BOOST_CHECK_EQUAL(read.duplicate(), original_read.duplicate());
  BOOST_CHECK_EQUAL(read.supplementary(), original_read.supplementary());
}

BOOST_AUTO_TEST_CASE( build_multiple_reads ) {
  auto header = SingleSamReader{"testdata/test_simple.bam"}.header();
  auto builder = SamBuilder{header};
  builder.set_name("foo").set_cigar("1M1I1M").set_bases("ACT").set_base_quals({1, 2, 3});
  builder.set_chromosome(1).set_alignment_start(142).set_mate_alignment_start(199);
  auto constructed_read = builder.build();

  BOOST_CHECK_EQUAL(constructed_read.name(), "foo");
  BOOST_CHECK_EQUAL(constructed_read.cigar().to_string(), "1M1I1M");
  BOOST_CHECK_EQUAL(constructed_read.bases().to_string(), "ACT");
  BOOST_CHECK_EQUAL(constructed_read.base_quals().to_string(), "1 2 3");
  BOOST_CHECK_EQUAL(constructed_read.chromosome(), 1u);
  BOOST_CHECK_EQUAL(constructed_read.alignment_start(), 142u);
  BOOST_CHECK_EQUAL(constructed_read.mate_alignment_start(), 199u);

  auto second_read = builder.set_chromosome(2).set_cigar("2M1I").set_bases("GGA").build();

  BOOST_CHECK_EQUAL(second_read.name(), "foo");
  BOOST_CHECK_EQUAL(second_read.cigar().to_string(), "2M1I");
  BOOST_CHECK_EQUAL(second_read.bases().to_string(), "GGA");
  BOOST_CHECK_EQUAL(second_read.base_quals().to_string(), "1 2 3");
  BOOST_CHECK_EQUAL(second_read.chromosome(), 2u);
  BOOST_CHECK_EQUAL(second_read.alignment_start(), 142u);
  BOOST_CHECK_EQUAL(second_read.mate_alignment_start(), 199u);
}

BOOST_AUTO_TEST_CASE( set_illegal_base_quals ) {
  auto header = SingleSamReader{"testdata/test_simple.bam"}.header();
  auto builder = SamBuilder{header};
  builder.set_name("foo").set_bases("ACT").set_cigar("3M");
  BOOST_CHECK_THROW(builder.set_base_quals({2, 255, 256}), invalid_argument); // 256 is too big for a quality score
}

BOOST_AUTO_TEST_CASE( set_illegal_tag ) {
  auto read = *(SingleSamReader{"testdata/test_simple.bam"}.begin());
  auto builder = SamBuilder{read};

  // Should throw an exception, since tag name is invalid
  BOOST_CHECK_THROW(builder.clear_tags().add_integer_tag("abc", 1).build(), invalid_argument);
  builder = SamBuilder{read};
  BOOST_CHECK_THROW(builder.clear_tags().add_double_tag("abc", 1.0).build(), invalid_argument);
  BOOST_CHECK_THROW(builder.clear_tags().add_char_tag("abc", 'A').build(), invalid_argument);
  BOOST_CHECK_THROW(builder.clear_tags().add_string_tag("abc", "A").build(), invalid_argument);

  // Integer tag has range of [-INT32, UINT32]
  BOOST_CHECK_THROW(builder.clear_tags().add_integer_tag("aa", LONG_MAX).build(), invalid_argument);
}

BOOST_AUTO_TEST_CASE( set_illegal_cigar ) {
  auto header = SingleSamReader{"testdata/test_simple.bam"}.header();
  auto builder = SamBuilder{header};
  builder.set_name("foo").set_bases("ACT").set_base_quals({1, 2, 3});

  // Should throw an exception, since the cigar is invalid
  BOOST_CHECK_THROW(builder.set_cigar("3F"), invalid_argument);// F is not a cigar operator
  BOOST_CHECK_THROW(builder.set_cigar("3"), invalid_argument); // missing operator
  BOOST_CHECK_THROW(builder.set_cigar("M3"), invalid_argument);// operator in front of the element length
  BOOST_CHECK_THROW(builder.set_cigar("3í"), invalid_argument);// operator is > 128 in ascii table
  BOOST_CHECK_THROW(builder.set_cigar("100M").build(), logic_error);// there are only 3 bases
}

BOOST_AUTO_TEST_CASE( set_mismatching_bases_and_quals ) {
  auto header = SingleSamReader{"testdata/test_simple.bam"}.header();
  auto builder = SamBuilder{header};
  builder.set_name("foo").set_cigar("3M").set_bases("ACT").set_base_quals({1, 2});

  // Should throw an exception, since the number of base qualities does not equal the number of bases
  BOOST_CHECK_THROW(builder.build(), logic_error);
}

BOOST_AUTO_TEST_CASE( set_mismatching_cigar_and_bases ) {
  auto header = SingleSamReader{"testdata/test_simple.bam"}.header();
  auto builder = SamBuilder{header};
  builder.set_name("foo").set_cigar("2M").set_bases("ACT").set_base_quals({1, 2, 3});

  // Should throw an exception, since the cigar does not match the sequence length
  BOOST_CHECK_THROW(builder.build(), logic_error);
}

BOOST_AUTO_TEST_CASE( missing_required_data_field ) {
  auto header = SingleSamReader{"testdata/test_simple.bam"}.header();
  auto builder = SamBuilder{header};
  builder.set_name("foo").set_bases("ACT").set_base_quals({1, 2, 3});

  // Should throw an exception, since we're missing a cigar
  BOOST_CHECK_THROW(builder.build(), logic_error);
}

BOOST_AUTO_TEST_CASE( build_without_validation ) {
  auto header = SingleSamReader{"testdata/test_simple.bam"}.header();
  auto builder = SamBuilder{header, false}; // disable validation
  builder.set_name("foo").set_cigar("2M").set_bases("ACT").set_base_quals({1, 2, 3});

  // Should not throw, since we've disabled validation
  auto read = builder.build();

  BOOST_CHECK_EQUAL(read.name(), "foo");
  BOOST_CHECK_EQUAL(read.cigar().to_string(), "2M");
  BOOST_CHECK_EQUAL(read.bases().to_string(), "ACT");
  BOOST_CHECK_EQUAL(read.base_quals().to_string(), "1 2 3");
}

BOOST_AUTO_TEST_CASE( builder_move_constructor ) {
  auto header = SingleSamReader{"testdata/test_simple.bam"}.header();
  auto builder1 = SamBuilder{header};
  const auto read1 = builder1.set_name("TEST_READ").set_cigar("1M").set_bases("A").set_base_quals(initializer_list<uint8_t>{24}).build();
  auto builder2 = std::move(builder1);
  const auto read2 = builder2.build();
  BOOST_CHECK_EQUAL(read1.name(), read2.name()); // the sheer fact that we can run this means the move constructor worked. Checking these just to make sure.
  BOOST_CHECK(read1.cigar()      == read2.cigar());
  BOOST_CHECK(read1.bases()      == read2.bases());
  BOOST_CHECK(read1.base_quals() == read2.base_quals());
  builder1 = std::move(builder2);     // check move assignment
  const auto read3 = builder1.build();
  BOOST_CHECK_EQUAL(read1.name(), read3.name()); // the sheer fact that we can run this means the move constructor worked. Checking these just to make sure.
  BOOST_CHECK(read1.cigar()      == read3.cigar());
  BOOST_CHECK(read1.bases()      == read3.bases());
  BOOST_CHECK(read1.base_quals() == read3.base_quals());
}

BOOST_AUTO_TEST_CASE( builder_starting_read_constructor ) {
  auto reader = SingleSamReader{"testdata/test_simple.bam"};
  auto original_read = *reader.begin();
  const auto header = reader.header();
  const auto builder = SamBuilder{header, original_read};
  original_read.set_unmapped();      // modifying the original read does not modify the reads created by the builder even during the build process.
  const auto read = builder.build(); // create a read (basically a copy of the original read)
  BOOST_CHECK(!read.unmapped());
}
