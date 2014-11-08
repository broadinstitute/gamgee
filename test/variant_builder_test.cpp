#include "variant_builder.h"
#include "variant_reader.h"
#include "missing.h"
#include "genotype.h"
#include "htslib/vcf.h"

#include <boost/test/unit_test.hpp>
#include <stdexcept>
#include <algorithm>

using namespace std;
using namespace gamgee;

/******************************************************
 * Core field tests (chromosome, alignment start, etc.)
 ******************************************************/

BOOST_AUTO_TEST_CASE( set_required_fields_only ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};

  // chromosome by string, one-base ref allele
  auto variant = builder.set_chromosome("20").set_alignment_start(5).set_ref_allele("A").build();
  BOOST_CHECK_EQUAL(variant.chromosome(), 1);
  BOOST_CHECK_EQUAL(variant.alignment_start(), 5);
  BOOST_CHECK_EQUAL(variant.ref(), "A");
  BOOST_CHECK_EQUAL(variant.alignment_stop(), 5);  // making sure rlen was set to length of ref allele here (1)

  BOOST_CHECK(missing(variant.qual()));
  BOOST_CHECK(missing(variant.id()));
  BOOST_CHECK(missing(variant.alt()));
  BOOST_CHECK_EQUAL(variant.filters().size(), 0);
  BOOST_CHECK(missing(variant.genotypes()));

  // chromosome by index, multi-base ref allele
  builder.clear();
  variant = builder.set_ref_allele("CA").set_chromosome(2).set_alignment_start(3).build();
  BOOST_CHECK_EQUAL(variant.ref(), "CA");
  BOOST_CHECK_EQUAL(variant.chromosome(), 2);
  BOOST_CHECK_EQUAL(variant.alignment_start(), 3);
  BOOST_CHECK_EQUAL(variant.alignment_stop(), 4);  // making sure rlen was set to length of ref allele here (2)

  BOOST_CHECK(missing(variant.qual()));
  BOOST_CHECK(missing(variant.id()));
  BOOST_CHECK(missing(variant.alt()));
  BOOST_CHECK_EQUAL(variant.filters().size(), 0);
  BOOST_CHECK(missing(variant.genotypes()));

  // long ref allele
  variant = builder.clear().set_chromosome(2).set_alignment_start(3).set_ref_allele("ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ").build();
  BOOST_CHECK_EQUAL(variant.ref(), "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ");
}

BOOST_AUTO_TEST_CASE( fail_to_set_required_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};

  // Failing to set one or more of chromosome, alignment start, or the ref allele is an error
  BOOST_CHECK_THROW(builder.build(), logic_error);
  builder.clear();
  BOOST_CHECK_THROW(builder.set_chromosome(0).build(), logic_error);
  builder.clear();
  BOOST_CHECK_THROW(builder.set_alignment_start(1).build(), logic_error);
  builder.clear();
  BOOST_CHECK_THROW(builder.set_ref_allele("A").build(), logic_error);
  builder.clear();
  BOOST_CHECK_THROW(builder.set_chromosome(0).set_alignment_start(1).build(), logic_error);
  builder.clear();
  BOOST_CHECK_THROW(builder.set_chromosome(0).set_ref_allele("A").build(), logic_error);
  builder.clear();
  BOOST_CHECK_THROW(builder.set_alignment_start(1).set_ref_allele("A").build(), logic_error);

  // Setting ref to a missing value should trigger an exception
  // (note that chromosome and alignment start can't be set to missing values, since
  //  they take unsigned arguments)
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  BOOST_CHECK_THROW(builder.set_ref_allele(""), invalid_argument);
  builder.set_ref_allele("A");
  BOOST_CHECK_THROW(builder.set_ref_allele("."), invalid_argument);
}

BOOST_AUTO_TEST_CASE( set_alignment_stop ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("ACG");

  // when not explicitly set, alignment stop should be based on ref allele length
  BOOST_CHECK_EQUAL(builder.build().alignment_stop(), 7);

  // alignment stop == alignment start
  auto variant = builder.set_alignment_stop(5).build();
  BOOST_CHECK_EQUAL(variant.alignment_stop(), 5);

  // alignment stop == alignment start + 1
  variant = builder.set_alignment_stop(6).build();
  BOOST_CHECK_EQUAL(variant.alignment_stop(), 6);

  // alignment stop > alignment start by more than 1
  variant = builder.set_alignment_stop(10).build();
  BOOST_CHECK_EQUAL(variant.alignment_stop(), 10);

  // alignment stop < alignment start is an error
  BOOST_CHECK_THROW(builder.set_alignment_stop(4).build(), logic_error);
  BOOST_CHECK_THROW(builder.set_alignment_stop(3).build(), logic_error);
  BOOST_CHECK_THROW(builder.set_alignment_stop(0).build(), logic_error);
}

BOOST_AUTO_TEST_CASE( set_qual ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Qual should be missing if not explicitly set
  BOOST_CHECK(missing(builder.build().qual()));

  // Set qual to positive value
  auto variant = builder.set_qual(1.5).build();
  BOOST_CHECK_EQUAL(variant.qual(), 1.5);

  // Set qual to 0.0 (should NOT be missing!)
  variant = builder.set_qual(0.0).build();
  BOOST_CHECK_EQUAL(variant.qual(), 0.0);
  BOOST_CHECK(! missing(variant.qual()));

  // Set qual to negative value
  variant = builder.set_qual(-1.5).build();
  BOOST_CHECK_EQUAL(variant.qual(), -1.5);

  // Setting qual to missing should work
  auto missing_val = 0.0f; bcf_float_set_missing(missing_val);
  variant = builder.set_qual(missing_val).build();
  BOOST_CHECK(missing(variant.qual()));
}

BOOST_AUTO_TEST_CASE( remove_core_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};

  // Clearing alignment stop via remove_alignment_stop() after setting it
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("AC").set_alignment_stop(10);
  BOOST_CHECK_EQUAL(builder.build().alignment_stop(), 10);
  builder.remove_alignment_stop();
  BOOST_CHECK_EQUAL(builder.build().alignment_stop(), 6);  // once cleared, should be based on ref allele length

  // Clearing qual via remove_qual() after setting it to non-missing
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A").set_qual(1.5);
  BOOST_CHECK_EQUAL(builder.build().qual(), 1.5);
  builder.remove_qual();
  BOOST_CHECK(missing(builder.build().qual()));
}

BOOST_AUTO_TEST_CASE( set_bad_contig ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};

  // Set non-existent contig by name
  builder.set_alignment_start(5).set_ref_allele("A");
  BOOST_CHECK_THROW(builder.set_chromosome("FOO").build(), invalid_argument);

  // Set non-existent contig by index
  BOOST_CHECK_THROW(builder.set_chromosome(500).build(), invalid_argument);
}

/************************************************************
 * Non-info shared field tests
 * (other than ref, which is a required field tested above)
 ************************************************************/

BOOST_AUTO_TEST_CASE( set_id ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // ID should be missing unless explicitly set
  BOOST_CHECK(missing(builder.build().id()));

  // setting ID to a non-missing value
  auto variant = builder.set_id("ABC").build();
  BOOST_CHECK_EQUAL(variant.id(), "ABC");

  // setting ID to a long value
  variant = builder.set_id("ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ").build();
  BOOST_CHECK_EQUAL(variant.id(), "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ");

  // Setting ID to a missing value without having previously set it to anything
  // (the case where it HAS been previously set is tested below in remove_non_info_shared_field)
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_id("");
  BOOST_CHECK(missing(builder.build().id()));
  builder.set_id(".");
  BOOST_CHECK(missing(builder.build().id()));
}

BOOST_AUTO_TEST_CASE( set_alt_alleles ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Alt alleles should be missing unless explicitly set
  BOOST_CHECK(missing(builder.build().alt()));

  // One alt allele, using the string setter
  auto variant = builder.set_alt_allele("C").build();
  auto alts = variant.alt();
  BOOST_CHECK_EQUAL(alts.size(), 1);
  BOOST_CHECK_EQUAL(alts[0], "C");

  // One alt allele, using the vector setter
  variant = builder.set_alt_alleles({"C"}).build();
  alts = variant.alt();
  BOOST_CHECK_EQUAL(alts.size(), 1);
  BOOST_CHECK_EQUAL(alts[0], "C");

  // Two alt alleles
  variant = builder.set_alt_alleles({"C", "ATG"}).build();
  alts = variant.alt();
  BOOST_CHECK_EQUAL(alts.size(), 2);
  BOOST_CHECK_EQUAL(alts[0], "C");
  BOOST_CHECK_EQUAL(alts[1], "ATG");

  // Three alt alleles
  variant = builder.set_alt_alleles({"C", "ATG", "T"}).build();
  alts = variant.alt();
  BOOST_CHECK_EQUAL(alts.size(), 3);
  BOOST_CHECK_EQUAL(alts[0], "C");
  BOOST_CHECK_EQUAL(alts[1], "ATG");
  BOOST_CHECK_EQUAL(alts[2], "T");

  // Long alt alleles
  variant = builder.set_alt_alleles({"CTGACTGACTGACTGACTGACTGACTGACTGATCGATCGATCGATCGCTAGCTAGCTCGATC", "TGCATGCTAGCTGATCGATCGGGGGATTCGAGGGCTTTAGGGCTA"}).build();
  alts = variant.alt();
  BOOST_CHECK_EQUAL(alts.size(), 2);
  BOOST_CHECK_EQUAL(alts[0], "CTGACTGACTGACTGACTGACTGACTGACTGATCGATCGATCGATCGCTAGCTAGCTCGATC");
  BOOST_CHECK_EQUAL(alts[1], "TGCATGCTAGCTGATCGATCGGGGGATTCGAGGGCTTTAGGGCTA");

  // Setting alt to a missing value without having previously set it to anything
  // (the case where it HAS been previously set is tested below in remove_non_info_shared_field)
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  variant = builder.set_alt_alleles({}).build();
  BOOST_CHECK(missing(variant.alt()));

  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  variant = builder.set_alt_allele("").build();
  BOOST_CHECK(missing(variant.alt()));

  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  variant = builder.set_alt_allele(".").build();
  BOOST_CHECK(missing(variant.alt()));

  // Setting a mix of non-missing and missing alts is an error
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  BOOST_CHECK_THROW(builder.set_alt_alleles({"C", "", "T"}), invalid_argument);
  BOOST_CHECK_THROW(builder.set_alt_alleles({"C", ".", "T"}), invalid_argument);
}

BOOST_AUTO_TEST_CASE( set_filters ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Filters should be missing unless explicitly set
  BOOST_CHECK(missing(builder.build().filters()));

  // Set one filter by name
  auto filters = builder.set_filters(vector<string>{"PASS"}).build().filters();
  BOOST_CHECK_EQUAL(filters.size(), 1);
  BOOST_CHECK_EQUAL(filters[0], "PASS");

  // Set one filter by index
  filters = builder.set_filters(vector<int32_t>{0}).build().filters();
  BOOST_CHECK_EQUAL(filters.size(), 1);
  BOOST_CHECK_EQUAL(filters[0], "PASS");

  // Set multiple filters by name
  filters = builder.set_filters(vector<string>{"LOW_QUAL", "MISSED"}).build().filters();
  BOOST_CHECK_EQUAL(filters.size(), 2);
  BOOST_CHECK_EQUAL(filters[0], "LOW_QUAL");
  BOOST_CHECK_EQUAL(filters[1], "MISSED");

  // Set multiple filters by index
  filters = builder.set_filters(vector<int32_t>{header.field_index("LOW_QUAL"), header.field_index("MISSED")}).build().filters();
  BOOST_CHECK_EQUAL(filters.size(), 2);
  BOOST_CHECK_EQUAL(filters[0], "LOW_QUAL");
  BOOST_CHECK_EQUAL(filters[1], "MISSED");

  // Set non-existent filter by name
  BOOST_CHECK_THROW(builder.set_filters(vector<string>{"FOO"}).build(), invalid_argument);

  // Set non-existent filter by index
  BOOST_CHECK_THROW(builder.set_filters(vector<int32_t>{1000}).build(), invalid_argument);

  // Set mix of valid and non-existent filters
  BOOST_CHECK_THROW(builder.set_filters(vector<string>{"PASS", "FOO", "LOW_QUAL"}), invalid_argument);
  BOOST_CHECK_THROW(builder.set_filters(vector<int32_t>{header.field_index("PASS"), 1000, header.field_index("LOW_QUAL")}), invalid_argument);

  // Set filters to missing without having previously set it
  // (the case where they HAVE been previously set is tested below in remove_non_info_shared_field)
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  BOOST_CHECK(missing(builder.set_filters(vector<int32_t>{}).build().filters()));

  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  BOOST_CHECK(missing(builder.set_filters(vector<string>{}).build().filters()));
}

BOOST_AUTO_TEST_CASE( remove_non_info_shared_field ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Remove pre-existing ID via remove_id()
  builder.set_id("ABC");
  builder.remove_id();
  BOOST_CHECK(missing(builder.build().id()));

  // Remove pre-existing ID by setting it to missing
  builder.set_id("ABC");
  builder.set_id("");
  BOOST_CHECK(missing(builder.build().id()));
  builder.set_id("ABC");
  builder.set_id(".");
  BOOST_CHECK(missing(builder.build().id()));

  // Remove pre-existing alts via remove_alt_alleles()
  builder.set_alt_alleles({"A", "CT"});
  builder.remove_alt_alleles();
  BOOST_CHECK(missing(builder.build().alt()));

  // Remove pre-existing alts by setting them to missing values
  builder.set_alt_alleles({"A", "CT"});
  builder.set_alt_alleles({});
  BOOST_CHECK(missing(builder.build().alt()));
  builder.set_alt_alleles({"A", "CT"});
  builder.set_alt_allele("");
  BOOST_CHECK(missing(builder.build().alt()));
  builder.set_alt_alleles({"A", "CT"});
  builder.set_alt_allele(".");
  BOOST_CHECK(missing(builder.build().alt()));

  // Remove pre-existing filters via remove_filters()
  builder.set_filters(vector<string>{"LOW_QUAL", "MISSED"});
  builder.remove_filters();
  BOOST_CHECK(missing(builder.build().filters()));

  // Remove pre-existing filters by setting them to missing
  builder.set_filters(vector<string>{"LOW_QUAL", "MISSED"});
  builder.set_filters(vector<string>{});
  BOOST_CHECK(missing(builder.build().filters()));
  builder.set_filters(vector<string>{"LOW_QUAL", "MISSED"});
  builder.set_filters(vector<int32_t>{});
  BOOST_CHECK(missing(builder.build().filters()));
}


/**********************************************
 * Helper functions for info/shared field tests
 **********************************************/

void check_integer_shared_field(const Variant& variant, const string& field, const vector<int32_t>& expected) {
  auto shared_field = variant.integer_shared_field(field);
  BOOST_CHECK_EQUAL(shared_field.size(), expected.size());
  for ( auto i = 0u; i < shared_field.size(); ++i ) {
    BOOST_CHECK_EQUAL(shared_field[i], expected[i]);
  }
}

void check_float_shared_field(const Variant& variant, const string& field, const vector<float>& expected) {
  auto shared_field = variant.float_shared_field(field);
  BOOST_CHECK_EQUAL(shared_field.size(), expected.size());
  for ( auto i = 0u; i < shared_field.size(); ++i ) {
    if ( bcf_float_is_missing(expected[i]) ) BOOST_CHECK(bcf_float_is_missing(shared_field[i]));
    else if ( bcf_float_is_vector_end(expected[i]) ) BOOST_CHECK(bcf_float_is_vector_end(shared_field[i]));
    else BOOST_CHECK_EQUAL(shared_field[i], expected[i]);
  }
}


/*******************************************
 * Info/Shared field tests
 *******************************************/

BOOST_AUTO_TEST_CASE( set_integer_shared_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  auto field_values = vector<int32_t>{};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Run each test below using both field indices and field names
  for ( bool use_field_index : { true, false } ) {

    // Single value, no vector
    auto variant = use_field_index ? builder.set_integer_shared_field(header.field_index("AN"), 5).build() : builder.set_integer_shared_field("AN", 5).build();
    check_integer_shared_field(variant, "AN", {5});

    // Single-element vector
    field_values = {5};
    variant = use_field_index ? builder.set_integer_shared_field(header.field_index("AN"), field_values).build() : builder.set_integer_shared_field("AN", field_values).build();
    check_integer_shared_field(variant, "AN", field_values);

    // Two-element vector
    field_values = {5, 10};
    variant = use_field_index ? builder.set_integer_shared_field(header.field_index("ZI"), field_values).build() : builder.set_integer_shared_field("ZI", field_values).build();
    check_integer_shared_field(variant, "ZI", field_values);

    // Four-element vector
    field_values = {5, 10, 15, 20};
    variant = use_field_index ? builder.set_integer_shared_field(header.field_index("ZI"), field_values).build() : builder.set_integer_shared_field("ZI", field_values).build();
    check_integer_shared_field(variant, "ZI", field_values);

    // Large vector
    field_values = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    variant = use_field_index ? builder.set_integer_shared_field(header.field_index("ZI"), field_values).build() : builder.set_integer_shared_field("ZI", field_values).build();
    check_integer_shared_field(variant, "ZI", field_values);

    // Multiple integer fields
    auto an_field_values = vector<int32_t>{3};
    auto zi_field_values = vector<int32_t>{5, 10, 15, 20};
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
    variant = use_field_index ? builder.set_integer_shared_field(header.field_index("AN"), an_field_values).set_integer_shared_field(header.field_index("ZI"), zi_field_values).build() :
                                builder.set_integer_shared_field("AN", an_field_values).set_integer_shared_field("ZI", zi_field_values).build();
    check_integer_shared_field(variant, "AN", an_field_values);
    check_integer_shared_field(variant, "ZI", zi_field_values);
  }

  // Set integer field to a missing value without having previously set it
  // (the case where you HAVE previously set it is tested in remove_info_shared_fields below)
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_shared_field("ZI", vector<int32_t>{});
  BOOST_CHECK(missing(builder.build().integer_shared_field("ZI")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_shared_field("ZI", missing_values::int32);
  BOOST_CHECK(missing(builder.build().integer_shared_field("ZI")));
}

BOOST_AUTO_TEST_CASE( set_float_shared_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  auto float_missing = 0.0f; bcf_float_set_missing(float_missing);
  auto field_values = vector<float>{};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A").set_alt_allele("G");

  // Run each test below using both field indices and field names
  for ( bool use_field_index : { true, false } ) {

    // Single value, no vector
    auto variant = use_field_index ? builder.set_float_shared_field(header.field_index("AF"), 5.5).build() : builder.set_float_shared_field("AF", 5.5).build();
    check_float_shared_field(variant, "AF", vector<float>{5.5});

    // Single-element vector
    field_values = {5.5};
    variant = use_field_index ? builder.set_float_shared_field(header.field_index("AF"), field_values).build() : builder.set_float_shared_field("AF", field_values).build();
    check_float_shared_field(variant, "AF", field_values);

    // Two-element vector
    field_values = {5.5, 10.5};
    variant = use_field_index ? builder.set_float_shared_field(header.field_index("ZF"), field_values).build() : builder.set_float_shared_field("ZF", field_values).build();
    check_float_shared_field(variant, "ZF", field_values);

    // Four-element vector
    field_values = {5.5, 10.5, 15.5, 20.5};
    variant = use_field_index ? builder.set_float_shared_field(header.field_index("ZF"), field_values).build() : builder.set_float_shared_field("ZF", field_values).build();
    check_float_shared_field(variant, "ZF", field_values);

    // Large vector
    field_values = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5};
    variant = use_field_index ? builder.set_float_shared_field(header.field_index("ZF"), field_values).build() : builder.set_float_shared_field("ZF", field_values).build();
    check_float_shared_field(variant, "ZF", field_values);

    // Multiple float fields
    auto af_field_values = vector<float>{3.5};
    auto zf_field_values = vector<float>{5.5, 10.5, 15.5, 20.5};
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
    variant = use_field_index ? builder.set_float_shared_field(header.field_index("AF"), af_field_values).set_float_shared_field(header.field_index("ZF"), zf_field_values).build() :
                                builder.set_float_shared_field("AF", af_field_values).set_float_shared_field("ZF", zf_field_values).build();
    check_float_shared_field(variant, "AF", af_field_values);
    check_float_shared_field(variant, "ZF", zf_field_values);
  }

  // Set float field to a missing value without having previously set it
  // (the case where you HAVE previously set it is tested in remove_info_shared_fields below)
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_float_shared_field("ZF", vector<float>{});
  BOOST_CHECK(missing(builder.build().float_shared_field("ZF")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_float_shared_field("ZF", float_missing);
  BOOST_CHECK(missing(builder.build().float_shared_field("ZF")));
}

BOOST_AUTO_TEST_CASE( set_string_shared_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Normal string value by field name
  auto variant = builder.set_string_shared_field("DESC", "helloworld").build();
  auto field = variant.string_shared_field("DESC");
  BOOST_CHECK_EQUAL(field.size(), 1);
  BOOST_CHECK_EQUAL(field[0], "helloworld");

  // Normal string value by field index
  variant = builder.set_string_shared_field(header.field_index("DESC"), "helloworld").build();
  field = variant.string_shared_field("DESC");
  BOOST_CHECK_EQUAL(field.size(), 1);
  BOOST_CHECK_EQUAL(field[0], "helloworld");

  // Long string value
  variant = builder.set_string_shared_field("DESC", "very long string that definitely absolutely positively will require tons and tons of memory").build();
  field = variant.string_shared_field("DESC");
  BOOST_CHECK_EQUAL(field.size(), 1);
  BOOST_CHECK_EQUAL(field[0], "very long string that definitely absolutely positively will require tons and tons of memory");

  // Multiple string fields
  variant = builder.set_string_shared_field("DESC", "helloworld").set_string_shared_field("ZS", "bjarne").build();
  auto desc_field = variant.string_shared_field("DESC");
  auto zs_field = variant.string_shared_field("ZS");
  BOOST_CHECK_EQUAL(desc_field.size(), 1);
  BOOST_CHECK_EQUAL(desc_field[0], "helloworld");
  BOOST_CHECK_EQUAL(zs_field.size(), 1);
  BOOST_CHECK_EQUAL(zs_field[0], "bjarne");

  // Set string field to a missing value without having previously set it
  // (the case where you HAVE previously set it is tested in remove_info_shared_fields below)
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_string_shared_field("ZS", "");
  BOOST_CHECK(missing(builder.build().string_shared_field("ZS")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_string_shared_field("ZS", ".");
  BOOST_CHECK(missing(builder.build().string_shared_field("ZS")));
}

BOOST_AUTO_TEST_CASE( set_boolean_shared_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};

  // Set boolean field by field name
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  auto variant = builder.set_boolean_shared_field("VALIDATED").build();
  BOOST_CHECK(variant.boolean_shared_field("VALIDATED"));

  // Set boolean field by field index
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  variant = builder.set_boolean_shared_field(header.field_index("VALIDATED")).build();
  BOOST_CHECK(variant.boolean_shared_field("VALIDATED"));

  // Set multiple boolean fields
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  variant = builder.set_boolean_shared_field("VALIDATED").set_boolean_shared_field("ZFLG").build();
  BOOST_CHECK(variant.boolean_shared_field("VALIDATED"));
  BOOST_CHECK(variant.boolean_shared_field("ZFLG"));
}

BOOST_AUTO_TEST_CASE( set_nonexistent_shared_field ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Non-existent field by name
  BOOST_CHECK_THROW(builder.set_integer_shared_field("FOOBAR", 1), invalid_argument);
  BOOST_CHECK_THROW(builder.set_float_shared_field("FOOBAR", 1.0), invalid_argument);
  BOOST_CHECK_THROW(builder.set_string_shared_field("FOOBAR", "hello"), invalid_argument);
  BOOST_CHECK_THROW(builder.set_boolean_shared_field("FOOBAR"), invalid_argument);

  // Non-existent field by index
  BOOST_CHECK_THROW(builder.set_integer_shared_field(1000, 1), invalid_argument);
  BOOST_CHECK_THROW(builder.set_float_shared_field(1000, 1.0), invalid_argument);
  BOOST_CHECK_THROW(builder.set_string_shared_field(1000, "hello"), invalid_argument);
  BOOST_CHECK_THROW(builder.set_boolean_shared_field(1000), invalid_argument);
}

BOOST_AUTO_TEST_CASE( set_shared_field_wrong_type ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Type mismatch by name
  BOOST_CHECK_THROW(builder.set_integer_shared_field("AF", 1), invalid_argument);
  BOOST_CHECK_THROW(builder.set_float_shared_field("AN", 1.5), invalid_argument);
  BOOST_CHECK_THROW(builder.set_string_shared_field("AF", "hello"), invalid_argument);
  BOOST_CHECK_THROW(builder.set_boolean_shared_field("AF"), invalid_argument);

  // Type mismatch by index
  BOOST_CHECK_THROW(builder.set_integer_shared_field(header.field_index("AF"), 1), invalid_argument);
  BOOST_CHECK_THROW(builder.set_float_shared_field(header.field_index("AN"), 1.5), invalid_argument);
  BOOST_CHECK_THROW(builder.set_string_shared_field(header.field_index("AF"), "hello"), invalid_argument);
  BOOST_CHECK_THROW(builder.set_boolean_shared_field(header.field_index("AF")), invalid_argument);
}

BOOST_AUTO_TEST_CASE( remove_integer_info_shared_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};

  // Remove integer shared field via remove_shared_field()
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_shared_field("ZI", vector<int32_t>{1, 2, 3});
  builder.remove_shared_field("ZI");
  BOOST_CHECK(missing(builder.build().integer_shared_field("ZI")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_shared_field("ZI", vector<int32_t>{1, 2, 3});
  builder.remove_shared_field(header.field_index("ZI"));
  BOOST_CHECK(missing(builder.build().integer_shared_field("ZI")));

  // Remove multiple integer shared fields via remove_shared_fields()
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_shared_field("AN", 1);
  builder.set_integer_shared_field("ZI", vector<int32_t>{1, 2, 3});
  builder.remove_shared_fields(vector<string>{"AN", "ZI"});
  auto variant = builder.build();
  BOOST_CHECK(missing(variant.integer_shared_field("AN")));
  BOOST_CHECK(missing(variant.integer_shared_field("ZI")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_shared_field("AN", 1);
  builder.set_integer_shared_field("ZI", vector<int32_t>{1, 2, 3});
  builder.remove_shared_fields(vector<uint32_t>{uint32_t(header.field_index("AN")), uint32_t(header.field_index("ZI"))});
  variant = builder.build();
  BOOST_CHECK(missing(variant.integer_shared_field("AN")));
  BOOST_CHECK(missing(variant.integer_shared_field("ZI")));

  // Remove integer shared field by setting it to a missing value
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_shared_field("ZI", vector<int32_t>{1, 2, 3});
  builder.set_integer_shared_field("ZI", vector<int32_t>{});
  BOOST_CHECK(missing(builder.build().integer_shared_field("ZI")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_shared_field("ZI", vector<int32_t>{1, 2, 3});
  builder.set_integer_shared_field("ZI", missing_values::int32);
  BOOST_CHECK(missing(builder.build().integer_shared_field("ZI")));
}

BOOST_AUTO_TEST_CASE( remove_float_info_shared_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  auto float_missing = 0.0f; bcf_float_set_missing(float_missing);

  // Remove float shared field via remove_shared_field()
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_float_shared_field("ZF", vector<float>{1.0, 2.0, 3.0});
  builder.remove_shared_field("ZF");
  BOOST_CHECK(missing(builder.build().float_shared_field("ZF")));

  // Remove multiple float shared fields via remove_shared_fields()
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_float_shared_field("AF", 1.0);
  builder.set_float_shared_field("ZF", vector<float>{1.0, 2.0, 3.0});
  builder.remove_shared_fields(vector<string>{"AF", "ZF"});
  auto variant = builder.build();
  BOOST_CHECK(missing(variant.float_shared_field("AF")));
  BOOST_CHECK(missing(variant.float_shared_field("ZF")));

  // Remove float shared field by setting it to a missing value
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_float_shared_field("ZF", vector<float>{1.0, 2.0, 3.0});
  builder.set_float_shared_field("ZF", vector<float>{});
  BOOST_CHECK(missing(builder.build().float_shared_field("ZF")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_float_shared_field("ZF", vector<float>{1.0, 2.0, 3.0});
  builder.set_float_shared_field("ZF", float_missing);
  BOOST_CHECK(missing(builder.build().float_shared_field("ZF")));
}

BOOST_AUTO_TEST_CASE( remove_string_info_shared_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};

  // Remove string shared field via remove_shared_field()
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_string_shared_field("ZS", "ABC");
  builder.remove_shared_field("ZS");
  BOOST_CHECK(missing(builder.build().string_shared_field("ZS")));

  // Remove multiple string shared fields via remove_shared_fields()
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_string_shared_field("DESC", "ABC");
  builder.set_string_shared_field("ZS", "FOO");
  builder.remove_shared_fields(vector<string>{"DESC", "ZS"});
  auto variant = builder.build();
  BOOST_CHECK(missing(variant.string_shared_field("DESC")));
  BOOST_CHECK(missing(variant.string_shared_field("ZS")));

  // Remove string shared field by setting it to a missing value
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_string_shared_field("ZS", "ABC");
  builder.set_string_shared_field("ZS", "");
  BOOST_CHECK(missing(builder.build().string_shared_field("ZS")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_string_shared_field("ZS", "ABC");
  builder.set_string_shared_field("ZS", ".");
  BOOST_CHECK(missing(builder.build().string_shared_field("ZS")));
}

BOOST_AUTO_TEST_CASE( remove_boolean_info_shared_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};

  // Remove boolean shared field via remove_shared_field()
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_boolean_shared_field("ZFLG");
  builder.remove_shared_field("ZFLG");
  BOOST_CHECK(missing(builder.build().boolean_shared_field("ZFLG")));

  // Remove multiple boolean shared fields via remove_shared_fields()
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_boolean_shared_field("VALIDATED");
  builder.set_boolean_shared_field("ZFLG");
  builder.remove_shared_fields(vector<string>{"VALIDATED", "ZFLG"});
  auto variant = builder.build();
  BOOST_CHECK(missing(variant.boolean_shared_field("VALIDATED")));
  BOOST_CHECK(missing(variant.boolean_shared_field("ZFLG")));
}

BOOST_AUTO_TEST_CASE( remove_subset_of_shared_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  auto zi_values = vector<int32_t>{1, 2, 3};
  auto zf_values = vector<float>{3.0, 4.5, 5.0};

  builder.set_integer_shared_field("ZI", zi_values);
  builder.set_float_shared_field("AF", 1.0);
  builder.set_float_shared_field("ZF", zf_values);
  builder.set_string_shared_field("ZS", "FOO");
  builder.set_boolean_shared_field("ZFLG");

  // Remove only a subset of the total shared fields set
  builder.remove_shared_fields(vector<string>{"AF", "ZS"});
  auto variant = builder.build();

  // And make sure only those fields were removed
  BOOST_CHECK(missing(variant.float_shared_field("AF")));
  BOOST_CHECK(missing(variant.string_shared_field("ZS")));

  check_integer_shared_field(variant, "ZI", zi_values);
  check_float_shared_field(variant, "ZF", zf_values);
  BOOST_CHECK(variant.boolean_shared_field("ZFLG"));
}

BOOST_AUTO_TEST_CASE( set_shared_field_after_removal ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  auto zi_final_values = vector<int32_t>{4, 5, 6, 7};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Need to make sure that the "removed" flag gets cleared when we re-set a field after removal
  builder.set_integer_shared_field("ZI", {1, 2, 3});
  builder.remove_shared_field("ZI");
  builder.set_integer_shared_field("ZI", zi_final_values);

  auto variant = builder.build();
  BOOST_CHECK( ! missing(variant.integer_shared_field("ZI")) );
  check_integer_shared_field(variant, "ZI", zi_final_values);
}

BOOST_AUTO_TEST_CASE( test_shared_field_memory_pool_compaction ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  auto zi_values = vector<int32_t>{1, 2, 3};
  auto zf_values = vector<float>{3.0, 4.5, 5.0};

  builder.set_string_shared_field("ZS", "FOO");
  builder.set_integer_shared_field("ZI", zi_values);
  builder.set_float_shared_field("ZF", zf_values);

  // Repeatedly set ZS in order to fill the shared memory pool multiple times and trigger
  // several rounds of memory compaction (VariantBuilderSharedRegion::compact_shared_buffer()).
  //
  // (Note that this is not going to be a very common occurrence -- typically users will clear()
  //  their builders between records instead of just setting repeatedly without clearing)
  for ( auto i = 0u; i < 30000; ++i ) {
    builder.set_string_shared_field("ZS", "BAR");
  }
  builder.set_string_shared_field("ZS", "BINGO");

  auto variant = builder.build();

  // Make sure the memory compaction process hasn't corrupted field values!
  BOOST_CHECK( ! missing(variant.integer_shared_field("ZI")) );
  BOOST_CHECK( ! missing(variant.float_shared_field("ZF")) );
  BOOST_CHECK( ! missing(variant.string_shared_field("ZS")) );

  check_integer_shared_field(variant, "ZI", zi_values);
  check_float_shared_field(variant, "ZF", zf_values);
  BOOST_CHECK_EQUAL(variant.string_shared_field("ZS")[0], "BINGO");
}

/*********************************************
 * Helper functions for individual field tests
 *********************************************/

void check_integer_individual_field(const Variant& variant, const string& field, const vector<vector<int32_t>>& expected) {
  auto sample_index = 0u;
  for ( const auto& sample_values : variant.integer_individual_field(field) ) {
    BOOST_CHECK_EQUAL(sample_values.size(), expected[sample_index].size());

    for ( auto value_index = 0u; value_index < sample_values.size(); ++value_index ) {
      BOOST_CHECK_EQUAL(sample_values[value_index], expected[sample_index][value_index]);
    }
    ++sample_index;
  }
}

void check_float_individual_field(const Variant& variant, const string& field, const vector<vector<float>>& expected) {
  auto sample_index = 0u;
  for ( const auto& sample_values : variant.float_individual_field(field) ) {
    BOOST_CHECK_EQUAL(sample_values.size(), expected[sample_index].size());

    for ( auto value_index = 0u; value_index < sample_values.size(); ++value_index ) {
      if ( bcf_float_is_missing(expected[sample_index][value_index]) ) BOOST_CHECK(bcf_float_is_missing(sample_values[value_index]));
      else if ( bcf_float_is_vector_end(expected[sample_index][value_index]) ) BOOST_CHECK(bcf_float_is_vector_end(sample_values[value_index]));
      else BOOST_CHECK_EQUAL(expected[sample_index][value_index], sample_values[value_index]);
    }
    ++sample_index;
  }
}

void check_string_individual_field(const Variant& variant, const string& field, const vector<string>& expected) {
  auto expected_max_length = 0u;
  for_each(expected.begin(), expected.end(), [&expected_max_length] (const string& str) { if ( str.length() > expected_max_length ) expected_max_length = str.length(); });

  auto sample_index = 0u;
  for ( const auto& sample_values : variant.string_individual_field(field) ) {
    BOOST_CHECK_EQUAL(sample_values.size(), expected_max_length);
    BOOST_CHECK_EQUAL(sample_values[0], expected[sample_index].empty() ? "." : expected[sample_index]);
    ++sample_index;
  }
}

void check_genotype_field(const Variant& variant, const vector<vector<int32_t>>& expected) {
  auto sample_index = 0u;
  for ( const auto& genotype : variant.genotypes() ) {
    BOOST_CHECK_EQUAL(genotype.size(), expected[sample_index].size());

    for ( auto allele = 0u; allele < genotype.size(); ++allele ) {
      BOOST_CHECK_EQUAL(genotype[allele], expected[sample_index][allele]);
    }
    ++sample_index;
  }
}

void init_integer_individual_field_by_sample(VariantBuilder& builder, const string& field, const vector<string>& samples,
                                             const vector<vector<int32_t>>& values, bool use_indices) {
  auto header = builder.header();
  auto sample_idx = 0;
  for ( const auto& sample : samples ) {
    // Don't set values for samples with no data
    if ( ! values[sample_idx].empty() ) {
      use_indices ? builder.set_integer_individual_field(header.field_index(field), sample_idx, values[sample_idx]) :
                    builder.set_integer_individual_field(field, sample, values[sample_idx]);
    }
    ++sample_idx;
  }
}

void init_integer_individual_field_by_sample(VariantBuilder& builder, const string& field, const vector<string>& samples,
                                             const vector<int32_t>& values, bool use_indices) {
  auto header = builder.header();
  auto sample_idx = 0;
  for ( const auto& sample : samples ) {
    use_indices ? builder.set_integer_individual_field(header.field_index(field), sample_idx, values[sample_idx]) :
                  builder.set_integer_individual_field(field, sample, values[sample_idx]);
    ++sample_idx;
  }
}

void init_float_individual_field_by_sample(VariantBuilder& builder, const string& field, const vector<string>& samples,
                                           const vector<vector<float>>& values, bool use_indices) {
  auto header = builder.header();
  auto sample_idx = 0;
  for ( const auto& sample : samples ) {
    // Don't set values for samples with no data
    if ( ! values[sample_idx].empty() ) {
      use_indices ? builder.set_float_individual_field(header.field_index(field), sample_idx, values[sample_idx]) :
                    builder.set_float_individual_field(field, sample, values[sample_idx]);
    }
    ++sample_idx;
  }
}

void init_float_individual_field_by_sample(VariantBuilder& builder, const string& field, const vector<string>& samples,
                                           const vector<float>& values, bool use_indices) {
  auto header = builder.header();
  auto sample_idx = 0;
  for ( const auto& sample : samples ) {
    use_indices ? builder.set_float_individual_field(header.field_index(field), sample_idx, values[sample_idx]) :
                  builder.set_float_individual_field(field, sample, values[sample_idx]);
    ++sample_idx;
  }
}

void init_string_individual_field_by_sample(VariantBuilder& builder, const string& field, const vector<string>& samples,
                                            const vector<string>& values, bool use_indices) {
  auto header = builder.header();
  auto sample_idx = 0;
  for ( const auto& sample : samples ) {
    // Don't set values for samples with no data
    if ( ! values[sample_idx].empty() ) {
      use_indices ? builder.set_string_individual_field(header.field_index(field), sample_idx, values[sample_idx]) :
                    builder.set_string_individual_field(field, sample, values[sample_idx]);
    }
    ++sample_idx;
  }
}

void init_genotype_field_by_sample(VariantBuilder& builder, const vector<string>& samples,
                                   const vector<vector<int32_t>>& values, bool use_indices) {
  auto sample_idx = 0;
  for ( const auto& sample : samples ) {
    // Don't set values for samples with no data
    if ( ! values[sample_idx].empty() ) {
      use_indices ? builder.set_genotype(sample_idx, values[sample_idx]) :
                    builder.set_genotype(sample, values[sample_idx]);
    }
    ++sample_idx;
  }
}


/**********************************************
 * Individual fields: tests for setting in bulk
 * (ie., setting an entire field at once)
 **********************************************/

BOOST_AUTO_TEST_CASE( bulk_set_integer_individual_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Run each test below with both lvalues and rvalues
  for ( bool perform_move : { true, false } ) {

    // Run each test below using both field indices and field names
    for ( bool use_field_index : { true, false } ) {

      auto flat_vector = vector<int32_t>{};
      auto nested_vector = vector<vector<int32_t>>{};
      auto expected = vector<vector<int32_t>>{};

      // Flat vector with one value per sample
      flat_vector = {1, 2, 3};
      expected = {{1}, {2}, {3}};
      auto variant = use_field_index ? (perform_move ? builder.set_integer_individual_field(header.field_index("ZIFMT"), move(flat_vector)).build() : builder.set_integer_individual_field(header.field_index("ZIFMT"), flat_vector).build()) :
                                       (perform_move ? builder.set_integer_individual_field("ZIFMT", move(flat_vector)).build() : builder.set_integer_individual_field("ZIFMT", flat_vector).build());
      check_integer_individual_field(variant, "ZIFMT", expected);

      // Flat vector with multiple values per sample
      flat_vector = {1, 2, 3, 4, 5, 6, 7, 8, 9};
      expected = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
      variant = use_field_index ? (perform_move ? builder.set_integer_individual_field(header.field_index("ZIFMT"), move(flat_vector)).build() : builder.set_integer_individual_field(header.field_index("ZIFMT"), flat_vector).build()) :
                                  (perform_move ? builder.set_integer_individual_field("ZIFMT", move(flat_vector)).build() : builder.set_integer_individual_field("ZIFMT", flat_vector).build());
      check_integer_individual_field(variant, "ZIFMT", expected);

      // Uniform-length nested vector
      nested_vector = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
      expected = nested_vector;
      variant = use_field_index ? (perform_move ? builder.set_integer_individual_field(header.field_index("ZIFMT"), move(nested_vector)).build() : builder.set_integer_individual_field(header.field_index("ZIFMT"), nested_vector).build()) :
                                  (perform_move ? builder.set_integer_individual_field("ZIFMT", move(nested_vector)).build() : builder.set_integer_individual_field("ZIFMT", nested_vector).build());
      check_integer_individual_field(variant, "ZIFMT", expected);

      // Non-uniform-length nested vector without any empty sample values
      nested_vector = {{1, 2, 3}, {4}, {5, 6}};
      expected = {{1, 2, 3}, {4, bcf_int32_vector_end, bcf_int32_vector_end}, {5, 6, bcf_int32_vector_end}};
      variant = use_field_index ? (perform_move ? builder.set_integer_individual_field(header.field_index("ZIFMT"), move(nested_vector)).build() : builder.set_integer_individual_field(header.field_index("ZIFMT"), nested_vector).build()) :
                                  (perform_move ? builder.set_integer_individual_field("ZIFMT", move(nested_vector)).build() : builder.set_integer_individual_field("ZIFMT", nested_vector).build());
      check_integer_individual_field(variant, "ZIFMT", expected);

      // Non-uniform-length nested vector with empty sample value
      nested_vector = {{1, 2, 3}, {4}, {}};
      expected = {{1, 2, 3}, {4, bcf_int32_vector_end, bcf_int32_vector_end}, {bcf_int32_missing, bcf_int32_vector_end, bcf_int32_vector_end}};
      variant = use_field_index ? (perform_move ? builder.set_integer_individual_field(header.field_index("ZIFMT"), move(nested_vector)).build() : builder.set_integer_individual_field(header.field_index("ZIFMT"), nested_vector).build()) :
                                  (perform_move ? builder.set_integer_individual_field("ZIFMT", move(nested_vector)).build() : builder.set_integer_individual_field("ZIFMT", nested_vector).build());
      check_integer_individual_field(variant, "ZIFMT", expected);

      // Non-uniform-length nested vector with empty sample value, int16 values
      nested_vector = {{30001, 30002, 30003}, {30004}, {}};
      expected = {{30001, 30002, 30003}, {30004, bcf_int32_vector_end, bcf_int32_vector_end}, {bcf_int32_missing, bcf_int32_vector_end, bcf_int32_vector_end}};
      variant = use_field_index ? (perform_move ? builder.set_integer_individual_field(header.field_index("ZIFMT"), move(nested_vector)).build() : builder.set_integer_individual_field(header.field_index("ZIFMT"), nested_vector).build()) :
                                  (perform_move ? builder.set_integer_individual_field("ZIFMT", move(nested_vector)).build() : builder.set_integer_individual_field("ZIFMT", nested_vector).build());
      check_integer_individual_field(variant, "ZIFMT", expected);

      // Non-uniform-length nested vector with empty sample value, int32 values
      nested_vector = {{70001, 70002, 70003}, {70004}, {}};
      expected = {{70001, 70002, 70003}, {70004, bcf_int32_vector_end, bcf_int32_vector_end}, {bcf_int32_missing, bcf_int32_vector_end, bcf_int32_vector_end}};
      variant = use_field_index ? (perform_move ? builder.set_integer_individual_field(header.field_index("ZIFMT"), move(nested_vector)).build() : builder.set_integer_individual_field(header.field_index("ZIFMT"), nested_vector).build()) :
                                  (perform_move ? builder.set_integer_individual_field("ZIFMT", move(nested_vector)).build() : builder.set_integer_individual_field("ZIFMT", nested_vector).build());
      check_integer_individual_field(variant, "ZIFMT", expected);
    }
  }
}

BOOST_AUTO_TEST_CASE( bulk_set_float_individual_fields ) {
  auto float_vector_end = 0.0f; bcf_float_set_vector_end(float_vector_end);
  auto float_missing = 0.0f; bcf_float_set_missing(float_missing);
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Run each test below with both lvalues and rvalues
  for ( bool perform_move : { true, false } ) {

    // Run each test below using both field indices and field names
    for ( bool use_field_index : { true, false } ) {

      auto flat_vector = vector<float>{};
      auto nested_vector = vector<vector<float>>{};
      auto expected = vector<vector<float>>{};

      // Flat vector with one value per sample
      flat_vector = {1.0, 2.0, 3.0};
      expected = {{1.0}, {2.0}, {3.0}};
      auto variant = use_field_index ? (perform_move ? builder.set_float_individual_field(header.field_index("ZFFMT"), move(flat_vector)).build() : builder.set_float_individual_field(header.field_index("ZFFMT"), flat_vector).build()) :
                                       (perform_move ? builder.set_float_individual_field("ZFFMT", move(flat_vector)).build() : builder.set_float_individual_field("ZFFMT", flat_vector).build());
      check_float_individual_field(variant, "ZFFMT", expected);

      // Flat vector with multiple values per sample
      flat_vector = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
      expected = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
      variant = use_field_index ? (perform_move ? builder.set_float_individual_field(header.field_index("ZFFMT"), move(flat_vector)).build() : builder.set_float_individual_field(header.field_index("ZFFMT"), flat_vector).build()) :
                                  (perform_move ? builder.set_float_individual_field("ZFFMT", move(flat_vector)).build() : builder.set_float_individual_field("ZFFMT", flat_vector).build());
      check_float_individual_field(variant, "ZFFMT", expected);

      // Uniform-length nested vector
      nested_vector = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
      expected = nested_vector;
      variant = use_field_index ? (perform_move ? builder.set_float_individual_field(header.field_index("ZFFMT"), move(nested_vector)).build() : builder.set_float_individual_field(header.field_index("ZFFMT"), nested_vector).build()) :
                                  (perform_move ? builder.set_float_individual_field("ZFFMT", move(nested_vector)).build() : builder.set_float_individual_field("ZFFMT", nested_vector).build());
      check_float_individual_field(variant, "ZFFMT", expected);

      // Non-uniform-length nested vector without any empty sample values
      nested_vector = {{1.0, 2.0, 3.0}, {4.0}, {5.0, 6.0}};
      expected = {{1.0, 2.0, 3.0}, {4.0, float_vector_end, float_vector_end}, {5.0, 6.0, float_vector_end}};
      variant = use_field_index ? (perform_move ? builder.set_float_individual_field(header.field_index("ZFFMT"), move(nested_vector)).build() : builder.set_float_individual_field(header.field_index("ZFFMT"), nested_vector).build()) :
                                  (perform_move ? builder.set_float_individual_field("ZFFMT", move(nested_vector)).build() : builder.set_float_individual_field("ZFFMT", nested_vector).build());
      check_float_individual_field(variant, "ZFFMT", expected);

      // Non-uniform-length nested vector with empty sample value
      nested_vector = {{1.0, 2.0, 3.0}, {4.0}, {}};
      expected = {{1.0, 2.0, 3.0}, {4.0, float_vector_end, float_vector_end}, {float_missing, float_vector_end, float_vector_end}};
      variant = use_field_index ? (perform_move ? builder.set_float_individual_field(header.field_index("ZFFMT"), move(nested_vector)).build() : builder.set_float_individual_field(header.field_index("ZFFMT"), nested_vector).build()) :
                                  (perform_move ? builder.set_float_individual_field("ZFFMT", move(nested_vector)).build() : builder.set_float_individual_field("ZFFMT", nested_vector).build());
      check_float_individual_field(variant, "ZFFMT", expected);
    }
  }
}

BOOST_AUTO_TEST_CASE( bulk_set_string_individual_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Run each test below with both lvalues and rvalues
  for ( bool perform_move : { true, false } ) {

    // Run each test below using both field indices and field names
    for ( bool use_field_index : { true, false } ) {

      auto flat_vector = vector<string>{};
      auto expected = vector<string>{};

      // Equal-length values for each sample
      flat_vector = {"abc", "def", "ghi"};
      expected = flat_vector;   // make a copy in case we perform a move
      auto variant = use_field_index ? (perform_move ? builder.set_string_individual_field(header.field_index("ZSFMT"), move(flat_vector)).build() : builder.set_string_individual_field(header.field_index("ZSFMT"), flat_vector).build()) :
                                       (perform_move ? builder.set_string_individual_field("ZSFMT", move(flat_vector)).build() : builder.set_string_individual_field("ZSFMT", flat_vector).build());
      check_string_individual_field(variant, "ZSFMT", expected);

      // Variable-length values across samples, no empty values
      flat_vector = {"abc", "d", "gh"};
      expected = flat_vector;   // make a copy in case we perform a move
      variant = use_field_index ? (perform_move ? builder.set_string_individual_field(header.field_index("ZSFMT"), move(flat_vector)).build() : builder.set_string_individual_field(header.field_index("ZSFMT"), flat_vector).build()) :
                                  (perform_move ? builder.set_string_individual_field("ZSFMT", move(flat_vector)).build() : builder.set_string_individual_field("ZSFMT", flat_vector).build());
      check_string_individual_field(variant, "ZSFMT", expected);

      // Variable-length values across samples, with an empty value
      flat_vector = {"abc", "", "gh"};
      expected = flat_vector;   // make a copy in case we perform a move
      variant = use_field_index ? (perform_move ? builder.set_string_individual_field(header.field_index("ZSFMT"), move(flat_vector)).build() : builder.set_string_individual_field(header.field_index("ZSFMT"), flat_vector).build()) :
                                  (perform_move ? builder.set_string_individual_field("ZSFMT", move(flat_vector)).build() : builder.set_string_individual_field("ZSFMT", flat_vector).build());
      check_string_individual_field(variant, "ZSFMT", expected);
    }
  }
}

BOOST_AUTO_TEST_CASE( bulk_set_genotype_field ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A").set_alt_alleles({"C", "T"});

  // Run each test below with both lvalues and rvalues
  for ( bool perform_move : { true, false } ) {

    auto flat_vector = vector<int32_t>{};
    auto nested_vector = vector<vector<int32_t>>{};
    auto expected = vector<vector<int32_t>>{};

    // Flat vector, one allele per sample
    flat_vector = {0, 1, 0}; Genotype::encode_genotype(flat_vector);
    expected = {{0}, {1}, {0}};
    auto variant = perform_move ? builder.set_genotypes(move(flat_vector)).build() :
                                  builder.set_genotypes(flat_vector).build();
    check_genotype_field(variant, expected);

    // Flat vector, two alleles per sample
    flat_vector = {0, 1, 1, 0, 0, 2}; Genotype::encode_genotype(flat_vector);
    expected = {{0, 1}, {1, 0}, {0, 2}};
    variant = perform_move ? builder.set_genotypes(move(flat_vector)).build() :
                             builder.set_genotypes(flat_vector).build();
    check_genotype_field(variant, expected);

    // Flat vector, three alleles per sample
    flat_vector = {0, 1, 2, 1, 0, 2, 0, 0, 2}; Genotype::encode_genotype(flat_vector);
    expected = {{0, 1, 2}, {1, 0, 2}, {0, 0, 2}};
    variant = perform_move ? builder.set_genotypes(move(flat_vector)).build() :
                             builder.set_genotypes(flat_vector).build();
    check_genotype_field(variant, expected);

    // Flat vector, varying ploidy with manual padding
    flat_vector = {0, 1, 2, 1, bcf_int32_vector_end, bcf_int32_vector_end, 0, 1, bcf_int32_vector_end}; Genotype::encode_genotype(flat_vector);
    expected = {{0, 1, 2}, {1, bcf_int32_vector_end, bcf_int32_vector_end}, {0, 1, bcf_int32_vector_end}};
    variant = perform_move ? builder.set_genotypes(move(flat_vector)).build() :
                             builder.set_genotypes(flat_vector).build();
    check_genotype_field(variant, expected);

    // Flat vector, varying ploidy with manual padding and a missing sample
    flat_vector = {0, 1, 2, -1, bcf_int32_vector_end, bcf_int32_vector_end, 0, 1, bcf_int32_vector_end}; Genotype::encode_genotype(flat_vector);
    expected = {{0, 1, 2}, {bcf_int32_missing, bcf_int32_vector_end, bcf_int32_vector_end}, {0, 1, bcf_int32_vector_end}};
    variant = perform_move ? builder.set_genotypes(move(flat_vector)).build() :
              builder.set_genotypes(flat_vector).build();
    check_genotype_field(variant, expected);

    // Nested vector, one allele per sample
    nested_vector = {{0}, {1}, {0}}; Genotype::encode_genotypes(nested_vector);
    expected = {{0}, {1}, {0}};
    variant = perform_move ? builder.set_genotypes(move(nested_vector)).build() :
                             builder.set_genotypes(nested_vector).build();
    check_genotype_field(variant, expected);

    // Nested vector, two alleles per sample
    nested_vector = {{0, 1}, {1, 0}, {0, 2}}; Genotype::encode_genotypes(nested_vector);
    expected = {{0, 1}, {1, 0}, {0, 2}};
    variant = perform_move ? builder.set_genotypes(move(nested_vector)).build() :
                             builder.set_genotypes(nested_vector).build();
    check_genotype_field(variant, expected);

    // Nested vector, three alleles per sample
    nested_vector = {{0, 1, 2}, {1, 0, 2}, {0, 0, 2}}; Genotype::encode_genotypes(nested_vector);
    expected = {{0, 1, 2}, {1, 0, 2}, {0, 0, 2}};
    variant = perform_move ? builder.set_genotypes(move(nested_vector)).build() :
                             builder.set_genotypes(nested_vector).build();
    check_genotype_field(variant, expected);

    // Nested vector, varying ploidy
    nested_vector = {{0, 1, 2}, {1, 0}, {0}}; Genotype::encode_genotypes(nested_vector);
    expected = {{0, 1, 2}, {1, 0, bcf_int32_vector_end}, {0, bcf_int32_vector_end, bcf_int32_vector_end}};
    variant = perform_move ? builder.set_genotypes(move(nested_vector)).build() :
                             builder.set_genotypes(nested_vector).build();
    check_genotype_field(variant, expected);

    // Nested vector, varying ploidy with missing sample
    nested_vector = {{0, 1, 2}, {1, 0}, {}}; Genotype::encode_genotypes(nested_vector);
    expected = {{0, 1, 2}, {1, 0, bcf_int32_vector_end}, {bcf_int32_missing, bcf_int32_vector_end, bcf_int32_vector_end}};
    variant = perform_move ? builder.set_genotypes(move(nested_vector)).build() :
                             builder.set_genotypes(nested_vector).build();
    check_genotype_field(variant, expected);
  }
}

BOOST_AUTO_TEST_CASE( bulk_set_individual_fields_wrong_number_of_values ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  BOOST_CHECK_EQUAL(header.n_samples(), 3);  // Make sure samples haven't been added to the test file,
                                             // as the tests below will be invalidated if they have

  // int field: less than num_samples values in flat vector
  BOOST_CHECK_THROW(builder.set_integer_individual_field("ZIFMT", {1, 2}), invalid_argument);
  // int field: number of values not divisible by num_samples in flat vector
  BOOST_CHECK_THROW(builder.set_integer_individual_field("ZIFMT", {1, 2, 3, 4}), invalid_argument);
  BOOST_CHECK_THROW(builder.set_integer_individual_field("ZIFMT", {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}), invalid_argument);
  // int field: length of first dimension not equal to num_samples in nested vector
  BOOST_CHECK_THROW(builder.set_integer_individual_field("ZIFMT", { {1, 2, 3}, {4, 5} }), invalid_argument);
  BOOST_CHECK_THROW(builder.set_integer_individual_field("ZIFMT", { {1, 2, 3}, {4, 5}, {6}, {7} }), invalid_argument);
  // int field: length of first dimension not equal to num_samples in nested vector, but divisible by num_samples
  BOOST_CHECK_THROW(builder.set_integer_individual_field("ZIFMT", { {1, 2, 3}, {4, 5}, {6}, {7}, {8}, {9} }), invalid_argument);

  // float field: less than num_samples values in flat vector
  BOOST_CHECK_THROW(builder.set_float_individual_field("ZFFMT", vector<float>{1.0, 2.0}), invalid_argument);
  // float field: number of values not divisible by num_samples in flat vector
  BOOST_CHECK_THROW(builder.set_float_individual_field("ZFFMT", vector<float>{1.0, 2.0, 3.0, 4.0}), invalid_argument);
  BOOST_CHECK_THROW(builder.set_float_individual_field("ZFFMT", vector<float>{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}), invalid_argument);
  // float field: length of first dimension not equal to num_samples in nested vector
  BOOST_CHECK_THROW(builder.set_float_individual_field("ZFFMT", vector<vector<float>>{ {1.0, 2.0, 3.0}, {4.0, 5.0} }), invalid_argument);
  BOOST_CHECK_THROW(builder.set_float_individual_field("ZFFMT", vector<vector<float>>{ {1.0, 2.0, 3.0}, {4.0, 5.0}, {6.0}, {7.0} }), invalid_argument);
  // float field: length of first dimension not equal to num_samples in nested vector, but divisible by num_samples
  BOOST_CHECK_THROW(builder.set_float_individual_field("ZFFMT", vector<vector<float>>{ {1.0, 2.0, 3.0}, {4.0, 5.0}, {6.0}, {7.0}, {8.0}, {9.0} }), invalid_argument);

  // string field: number of values < num_samples
  BOOST_CHECK_THROW(builder.set_string_individual_field("ZSFMT", {"abc", "def"}), invalid_argument);;
  // string field: number of values > num_samples
  BOOST_CHECK_THROW(builder.set_string_individual_field("ZSFMT", {"abc", "d", "e", "f"}), invalid_argument);
  // string field: number of values > num_samples, but divisible by num_samples
  BOOST_CHECK_THROW(builder.set_string_individual_field("ZSFMT", {"abc", "d", "e", "f", "g", "h"}), invalid_argument);
}

BOOST_AUTO_TEST_CASE( bulk_set_individual_field_type_mismatch ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Accessing float field as an integer field
  BOOST_CHECK_THROW(builder.set_integer_individual_field("ZFFMT", {1, 2, 3}), invalid_argument);
  BOOST_CHECK_THROW(builder.set_integer_individual_field(header.field_index("ZFFMT"), {1, 2, 3}), invalid_argument);

  // Accessing integer field as a float field
  BOOST_CHECK_THROW(builder.set_float_individual_field("ZIFMT", {1.5f, 2.5f, 3.5f}), invalid_argument);
  BOOST_CHECK_THROW(builder.set_float_individual_field(header.field_index("ZIFMT"), {1.5f, 2.5f, 3.5f}), invalid_argument);

  // Accessing integer field as a string field
  BOOST_CHECK_THROW(builder.set_string_individual_field("ZIFMT", {"a", "b", "c"}), invalid_argument);
  BOOST_CHECK_THROW(builder.set_string_individual_field(header.field_index("ZIFMT"), {"a", "b", "c"}), invalid_argument);
}

BOOST_AUTO_TEST_CASE( bulk_set_individual_field_bad_field_id ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  BOOST_CHECK_THROW(builder.set_integer_individual_field("FOOBAR", {1, 2, 3}), invalid_argument);   // bad field name
  BOOST_CHECK_THROW(builder.set_integer_individual_field(500, {1, 2, 3}), invalid_argument);        // bad field index

  BOOST_CHECK_THROW(builder.set_float_individual_field("FOOBAR", {1.5f, 2.5f, 3.5f}), invalid_argument);  // bad field name
  BOOST_CHECK_THROW(builder.set_float_individual_field(500, {1.5f, 2.5f, 3.5f}), invalid_argument);       // bad field index

  BOOST_CHECK_THROW(builder.set_string_individual_field("FOOBAR", {"a", "b", "c"}), invalid_argument);   // bad field name
  BOOST_CHECK_THROW(builder.set_string_individual_field(500, {"a", "b", "c"}), invalid_argument);        // bad field index
}

BOOST_AUTO_TEST_CASE( bulk_set_individual_field_all_empty_values ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Single-dimensional empty vectors
  auto variant = builder.set_integer_individual_field("ZIFMT", vector<int32_t>{}).build();
  BOOST_CHECK(missing(variant.integer_individual_field("ZIFMT")));

  variant = builder.set_float_individual_field("ZFFMT", vector<float>{}).build();
  BOOST_CHECK(missing(variant.float_individual_field("ZFFMT")));

  variant = builder.set_string_individual_field("ZSFMT", vector<string>{}).build();
  BOOST_CHECK(missing(variant.string_individual_field("ZSFMT")));

  // Multi-dimensional empty vectors
  variant = builder.set_integer_individual_field("ZIFMT", vector<vector<int32_t>>{}).build();
  BOOST_CHECK(missing(variant.integer_individual_field("ZIFMT")));

  variant = builder.set_float_individual_field("ZFFMT", vector<vector<float>>{}).build();
  BOOST_CHECK(missing(variant.float_individual_field("ZFFMT")));

  // Multi-dimensional NON-empty vectors with ALL empty inner vectors
  variant = builder.set_integer_individual_field("ZIFMT", vector<vector<int32_t>>{ {}, {}, {} }).build();
  BOOST_CHECK(missing(variant.integer_individual_field("ZIFMT")));

  variant = builder.set_float_individual_field("ZFFMT", vector<vector<float>>{ {}, {}, {} }).build();
  BOOST_CHECK(missing(variant.float_individual_field("ZFFMT")));

  // String vector with ALL empty strings
  variant = builder.set_string_individual_field("ZSFMT", vector<string>{ "", "", "" }).build();
  BOOST_CHECK(missing(variant.string_individual_field("ZSFMT")));
}

BOOST_AUTO_TEST_CASE( remove_integer_individual_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};

  // Remove integer individual field via remove_individual_field()
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_individual_field("ZIFMT", {1, 2, 3});
  builder.remove_individual_field("ZIFMT");
  BOOST_CHECK(missing(builder.build().integer_individual_field("ZIFMT")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_individual_field("ZIFMT", {1, 2, 3});
  builder.remove_individual_field(header.field_index("ZIFMT"));
  BOOST_CHECK(missing(builder.build().integer_individual_field("ZIFMT")));

  // Remove multiple integer individual fields via remove_individual_fields()
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_individual_field("ZIFMT", {1, 2, 3});
  builder.set_integer_individual_field("GQ", {4, 5, 6});
  builder.remove_individual_fields(vector<string>{"ZIFMT", "GQ"});
  auto variant = builder.build();
  BOOST_CHECK(missing(variant.integer_individual_field("ZIFMT")));
  BOOST_CHECK(missing(variant.integer_individual_field("GQ")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_individual_field("ZIFMT", {1, 2, 3});
  builder.set_integer_individual_field("GQ", {4, 5, 6});
  builder.remove_individual_fields(vector<uint32_t>{uint32_t(header.field_index("ZIFMT")), uint32_t(header.field_index("GQ"))});
  variant = builder.build();
  BOOST_CHECK(missing(variant.integer_individual_field("ZIFMT")));
  BOOST_CHECK(missing(variant.integer_individual_field("GQ")));

  // Remove integer individual field by setting it to a missing value
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_individual_field("ZIFMT", {1, 2, 3});
  builder.set_integer_individual_field("ZIFMT", vector<int32_t>{});
  BOOST_CHECK(missing(builder.build().integer_individual_field("ZIFMT")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_individual_field("ZIFMT", {1, 2, 3});
  builder.set_integer_individual_field("ZIFMT", vector<vector<int32_t>>{});
  BOOST_CHECK(missing(builder.build().integer_individual_field("ZIFMT")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_individual_field("ZIFMT", {1, 2, 3});
  builder.set_integer_individual_field("ZIFMT", vector<vector<int32_t>>{ {}, {}, {} });
  BOOST_CHECK(missing(builder.build().integer_individual_field("ZIFMT")));
}

BOOST_AUTO_TEST_CASE( remove_float_individual_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};

  // Remove float individual field via remove_individual_field()
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_float_individual_field("ZFFMT", vector<float>{1.0, 2.0, 3.0});
  builder.remove_individual_field("ZFFMT");
  BOOST_CHECK(missing(builder.build().float_individual_field("ZFFMT")));

  // Remove multiple float individual fields via remove_individual_fields()
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_float_individual_field("ZFFMT", vector<float>{1.0, 2.0, 3.0});
  builder.set_float_individual_field("AF", vector<float>{1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
  builder.remove_individual_fields(vector<string>{"ZFFMT", "AF"});
  auto variant = builder.build();
  BOOST_CHECK(missing(variant.float_individual_field("ZFFMT")));
  BOOST_CHECK(missing(variant.float_individual_field("AF")));

  // Remove float individual field by setting it to a missing value
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_float_individual_field("ZFFMT", vector<float>{1.0, 2.0, 3.0});
  builder.set_float_individual_field("ZFFMT", vector<float>{});
  BOOST_CHECK(missing(builder.build().float_individual_field("ZFFMT")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_float_individual_field("ZFFMT", vector<float>{1.0, 2.0, 3.0});
  builder.set_float_individual_field("ZFFMT", vector<vector<float>>{});
  BOOST_CHECK(missing(builder.build().float_individual_field("ZFFMT")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_float_individual_field("ZFFMT", vector<float>{1.0, 2.0, 3.0});
  builder.set_float_individual_field("ZFFMT", vector<vector<float>>{ {}, {}, {} });
  BOOST_CHECK(missing(builder.build().float_individual_field("ZFFMT")));
}

BOOST_AUTO_TEST_CASE( remove_string_individual_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};

  // Remove string individual field via remove_individual_field()
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_string_individual_field("ZSFMT", vector<string>{"A", "B", "C"});
  builder.remove_individual_field("ZSFMT");
  BOOST_CHECK(missing(builder.build().string_individual_field("ZSFMT")));

  // Remove multiple string individual fields via remove_individual_fields()
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_string_individual_field("ZSFMT", vector<string>{"A", "B", "C"});
  builder.set_string_individual_field("AS", vector<string>{"D", "E", "F"});
  builder.remove_individual_fields(vector<string>{"ZSFMT", "AS"});
  auto variant = builder.build();
  BOOST_CHECK(missing(variant.string_individual_field("ZSFMT")));
  BOOST_CHECK(missing(variant.string_individual_field("AS")));

  // Remove string individual field by setting it to a missing value
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_string_individual_field("ZSFMT", vector<string>{"A", "B", "C"});
  builder.set_string_individual_field("ZSFMT", vector<string>{});
  BOOST_CHECK(missing(builder.build().string_individual_field("ZSFMT")));
  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_string_individual_field("ZSFMT", vector<string>{"A", "B", "C"});
  builder.set_string_individual_field("ZSFMT", vector<string>{ "", "", "" });
  BOOST_CHECK(missing(builder.build().string_individual_field("ZSFMT")));
}

BOOST_AUTO_TEST_CASE( remove_subset_of_individual_fields ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  auto encoded_gts = vector<int32_t>{0, 1, 1, 0, 0, 2}; Genotype::encode_genotype(encoded_gts);
  auto gt_values = vector<vector<int32_t>>{{0, 1}, {1, 0}, {0, 2}};
  auto zifmt_values = vector<vector<int32_t>>{{1}, {2}, {3}};
  auto zffmt_values = vector<vector<float>>{{1.0}, {2.0}, {3.0}};

  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A").set_alt_alleles({"T", "C"});
  builder.set_genotypes(encoded_gts);
  builder.set_integer_individual_field("ZIFMT", zifmt_values);
  builder.set_integer_individual_field("GQ", {4, 5, 6});
  builder.set_float_individual_field("ZFFMT", zffmt_values);
  builder.set_string_individual_field("ZSFMT", vector<string>{"A", "B", "C"});

  // Remove only a subset of the total individual fields set
  builder.remove_individual_fields(vector<string>{"GQ", "ZSFMT"});
  auto variant = builder.build();

  // And make sure only those individual fields were removed
  BOOST_CHECK(missing(variant.integer_individual_field("GQ")));
  BOOST_CHECK(missing(variant.string_individual_field("ZSFMT")));

  check_genotype_field(variant, gt_values);
  check_integer_individual_field(variant, "ZIFMT", zifmt_values);
  check_float_individual_field(variant, "ZFFMT", zffmt_values);
}

BOOST_AUTO_TEST_CASE( set_individual_field_after_removal ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  auto zifmt_final_values = vector<vector<int32_t>>{{4, 5}, {6, 7}, {8, 9}};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Need to make sure that the "removed" flag gets cleared when we re-set a field after removal
  builder.set_integer_individual_field("ZIFMT", {1, 2, 3});
  builder.remove_individual_field("ZIFMT");
  builder.set_integer_individual_field("ZIFMT", zifmt_final_values);

  auto variant = builder.build();
  BOOST_CHECK( ! missing(variant.integer_individual_field("ZIFMT")) );
  check_integer_individual_field(variant, "ZIFMT", zifmt_final_values);
}

/*************************************************
 * Individual fields: tests for setting by sample
 *************************************************/

BOOST_AUTO_TEST_CASE( set_integer_individual_fields_by_sample ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  auto sample_names = vector<string>{ "NA12878", "NA12891", "NA12892" };
  auto expected = vector<vector<int32_t>>{};

  // Run each test below using both field/sample indices and field/sample names
  for ( bool use_indices : { true, false } ) {

    // Set single value per sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
    auto per_sample_single_values = vector<int32_t>{ 4, 5, 6 };
    expected = { {4}, {5}, {6} };
    init_integer_individual_field_by_sample(builder, "ZIFMT", sample_names, per_sample_single_values, use_indices);
    auto variant = builder.build();
    check_integer_individual_field(variant, "ZIFMT", expected);

    // Set multiple values per sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
    auto per_sample_multi_values = vector<vector<int32_t>>{ {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };
    expected = per_sample_multi_values;
    init_integer_individual_field_by_sample(builder, "ZIFMT", sample_names, per_sample_multi_values, use_indices);
    variant = builder.build();
    check_integer_individual_field(variant, "ZIFMT", expected);

    // Set unequal number of values per sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
    auto per_sample_ragged_values = vector<vector<int32_t>>{ {1, 2, 3}, {4}, {5, 6} };
    expected = { {1, 2, 3}, {4, bcf_int32_vector_end, bcf_int32_vector_end}, {5, 6, bcf_int32_vector_end} };
    init_integer_individual_field_by_sample(builder, "ZIFMT", sample_names, per_sample_ragged_values, use_indices);
    variant = builder.build();
    check_integer_individual_field(variant, "ZIFMT", expected);

    // Set unequal number of values per sample, with a missing sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
    auto per_sample_ragged_values_with_missing_sample = vector<vector<int32_t>>{ {1, 2, 3}, {}, {4, 5} };
    expected = { {1, 2, 3}, {bcf_int32_missing, bcf_int32_vector_end, bcf_int32_vector_end}, {4, 5, bcf_int32_vector_end} };
    init_integer_individual_field_by_sample(builder, "ZIFMT", sample_names, per_sample_ragged_values_with_missing_sample, use_indices);
    variant = builder.build();
    check_integer_individual_field(variant, "ZIFMT", expected);
  }
}

BOOST_AUTO_TEST_CASE( set_float_individual_fields_by_sample ) {
  auto float_missing = 0.0f;  bcf_float_set_missing(float_missing);
  auto float_vector_end = 0.0f;  bcf_float_set_vector_end(float_vector_end);
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  auto sample_names = vector<string>{ "NA12878", "NA12891", "NA12892" };
  auto expected = vector<vector<float>>{};

  // Run each test below using both field/sample indices and field/sample names
  for ( bool use_indices : { true, false } ) {

    // Set single value per sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
    auto per_sample_single_values = vector<float>{ 4.5, 5.5, 6.5 };
    expected = { {4.5}, {5.5}, {6.5} };
    init_float_individual_field_by_sample(builder, "ZFFMT", sample_names, per_sample_single_values, use_indices);
    auto variant = builder.build();
    check_float_individual_field(variant, "ZFFMT", expected);

    // Set multiple values per sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
    auto per_sample_multi_values = vector<vector<float>>{ {1.5, 2.5, 3.5}, {4.5, 5.5, 6.5}, {7.5, 8.5, 9.5} };
    expected = per_sample_multi_values;
    init_float_individual_field_by_sample(builder, "ZFFMT", sample_names, per_sample_multi_values, use_indices);
    variant = builder.build();
    check_float_individual_field(variant, "ZFFMT", expected);

    // Set unequal number of values per sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
    auto per_sample_ragged_values = vector<vector<float>>{ {1.5, 2.5, 3.5}, {4.5}, {5.5, 6.5} };
    expected = vector<vector<float>>{ {1.5, 2.5, 3.5}, {4.5, float_vector_end, float_vector_end}, {5.5, 6.5, float_vector_end} };
    init_float_individual_field_by_sample(builder, "ZFFMT", sample_names, per_sample_ragged_values, use_indices);
    variant = builder.build();
    check_float_individual_field(variant, "ZFFMT", expected);

    // Set unequal number of values per sample, with a missing sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
    auto per_sample_ragged_values_with_missing_sample = vector<vector<float>>{ {1.5, 2.5, 3.5}, {}, {4.5, 5.5} };
    expected = { {1.5, 2.5, 3.5}, {float_missing, float_vector_end, float_vector_end}, {4.5, 5.5, float_vector_end} };
    init_float_individual_field_by_sample(builder, "ZFFMT", sample_names, per_sample_ragged_values_with_missing_sample, use_indices);
    variant = builder.build();
    check_float_individual_field(variant, "ZFFMT", expected);
  }
}

BOOST_AUTO_TEST_CASE( set_string_individual_fields_by_sample ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  auto sample_names = vector<string>{ "NA12878", "NA12891", "NA12892" };
  auto expected = vector<string>{};

  // Run each test below using both field/sample indices and field/sample names
  for ( bool use_indices : { true, false } ) {

    // Set equal-length values per sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
    auto per_sample_equal_values = vector<string>{ "abc", "def", "ghi" };
    expected = per_sample_equal_values;
    init_string_individual_field_by_sample(builder, "ZSFMT", sample_names, per_sample_equal_values, use_indices);
    auto variant = builder.build();
    check_string_individual_field(variant, "ZSFMT", expected);

    // Set unequal-length values per sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
    auto per_sample_ragged_values = vector<string>{ "abc", "d", "ef" };
    expected = per_sample_ragged_values;
    init_string_individual_field_by_sample(builder, "ZSFMT", sample_names, per_sample_ragged_values, use_indices);
    variant = builder.build();
    check_string_individual_field(variant, "ZSFMT", expected);

    // Set unequal number of values per sample, with a missing sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
    auto per_sample_ragged_values_with_missing_sample = vector<string>{ "abc", "", "de" };
    expected = per_sample_ragged_values_with_missing_sample;
    init_string_individual_field_by_sample(builder, "ZSFMT", sample_names, per_sample_ragged_values_with_missing_sample, use_indices);
    variant = builder.build();
    check_string_individual_field(variant, "ZSFMT", expected);
  }
}

BOOST_AUTO_TEST_CASE( set_genotype_field_by_sample ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  auto sample_names = vector<string>{ "NA12878", "NA12891", "NA12892" };
  auto genotype_values = vector<vector<int32_t>>{};
  auto expected = vector<vector<int32_t>>{};

  // Run each test below using both field/sample indices and field/sample names
  for ( bool use_indices : { true, false } ) {

    // One allele per sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A").set_alt_alleles({"C", "T"});
    genotype_values = { {1}, {2}, {0} }; Genotype::encode_genotypes(genotype_values);
    expected = { {1}, {2}, {0} };
    init_genotype_field_by_sample(builder, sample_names, genotype_values, use_indices);
    auto variant = builder.build();
    check_genotype_field(variant, expected);

    // Two alleles per sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A").set_alt_alleles({"C", "T"});
    genotype_values = { {1, 2}, {0, 1}, {1, 1} }; Genotype::encode_genotypes(genotype_values);
    expected = { {1, 2}, {0, 1}, {1, 1} };
    init_genotype_field_by_sample(builder, sample_names, genotype_values, use_indices);
    variant = builder.build();
    check_genotype_field(variant, expected);

    // Three alleles per sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A").set_alt_alleles({"C", "T"});
    genotype_values = { {1, 0, 0}, {2, 1, 0}, {0, 2, 0} }; Genotype::encode_genotypes(genotype_values);
    expected = { {1, 0, 0}, {2, 1, 0}, {0, 2, 0} };
    init_genotype_field_by_sample(builder, sample_names, genotype_values, use_indices);
    variant = builder.build();
    check_genotype_field(variant, expected);

    // Varying ploidy, no missing samples
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A").set_alt_alleles({"C", "T"});
    genotype_values = { {0, 1, 2}, {1, 0}, {0} }; Genotype::encode_genotypes(genotype_values);
    expected = { {0, 1, 2}, {1, 0, bcf_int32_vector_end}, {0, bcf_int32_vector_end, bcf_int32_vector_end} };
    init_genotype_field_by_sample(builder, sample_names, genotype_values, use_indices);
    variant = builder.build();
    check_genotype_field(variant, expected);

    // Varying ploidy, with missing sample
    builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A").set_alt_alleles({"C", "T"});
    genotype_values = { {0, 1, 2}, {1, 0}, {} }; Genotype::encode_genotypes(genotype_values);
    expected = { {0, 1, 2}, {1, 0, bcf_int32_vector_end}, {bcf_int32_missing, bcf_int32_vector_end, bcf_int32_vector_end} };
    init_genotype_field_by_sample(builder, sample_names, genotype_values, use_indices);
    variant = builder.build();
    check_genotype_field(variant, expected);
  }
}

BOOST_AUTO_TEST_CASE( set_individual_field_by_sample_all_empty_values ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  auto sample_names = vector<string>{ "NA12878", "NA12891", "NA12892" };
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A").set_alt_alleles({"T", "C"});

  // Setting a field to ALL empty values by sample should have no effect
  // (field should not even be present in the final variant)
  for ( const auto& sample : sample_names ) {
    builder.set_integer_individual_field("ZIFMT", sample, vector<int32_t>{});
    builder.set_float_individual_field("ZFFMT", sample, vector<float>{});
    builder.set_string_individual_field("ZSFMT", sample, "");
    builder.set_genotype(sample, vector<int32_t>{});
  }
  auto variant = builder.build();

  BOOST_CHECK(missing(variant.integer_individual_field("ZIFMT")));
  BOOST_CHECK(missing(variant.float_individual_field("ZFFMT")));
  BOOST_CHECK(missing(variant.string_individual_field("ZSFMT")));
  BOOST_CHECK(missing(variant.genotypes()));
}

BOOST_AUTO_TEST_CASE( set_individual_field_by_sample_bad_sample_id ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  BOOST_CHECK_THROW(builder.set_integer_individual_field("ZIFMT", "FOOSAMPLE", {1, 2}), invalid_argument);             // bad sample name
  BOOST_CHECK_THROW(builder.set_integer_individual_field(header.field_index("ZIFMT"), 50, {1, 2}), invalid_argument);  // bad sample index

  BOOST_CHECK_THROW(builder.set_float_individual_field("ZFFMT", "FOOSAMPLE", {1.5f, 2.5f}), invalid_argument);             // bad sample name
  BOOST_CHECK_THROW(builder.set_float_individual_field(header.field_index("ZFFMT"), 50, {1.5f, 2.5f}), invalid_argument);  // bad sample index

  BOOST_CHECK_THROW(builder.set_string_individual_field("ZSFMT", "FOOSAMPLE", "abc"), invalid_argument);             // bad sample name
  BOOST_CHECK_THROW(builder.set_string_individual_field(header.field_index("ZSFMT"), 50, "abc"), invalid_argument);  // bad sample index
}

BOOST_AUTO_TEST_CASE( set_individual_field_by_sample_bad_field_id ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  auto sample_names = vector<string>{ "NA12878", "NA12891", "NA12892" };
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  BOOST_CHECK_THROW(builder.set_integer_individual_field("FOOBAR", sample_names[0], {1, 2}), invalid_argument);        // bad field name
  BOOST_CHECK_THROW(builder.set_integer_individual_field(500, 0, {1, 2}), invalid_argument);                           // bad field index

  BOOST_CHECK_THROW(builder.set_float_individual_field("FOOBAR", sample_names[0], {1.5f, 2.5f}), invalid_argument);        // bad field name
  BOOST_CHECK_THROW(builder.set_float_individual_field(500, 0, {1.5f, 2.5f}), invalid_argument);                           // bad field index

  BOOST_CHECK_THROW(builder.set_string_individual_field("FOOBAR", sample_names[0], "abc"), invalid_argument);        // bad field name
  BOOST_CHECK_THROW(builder.set_string_individual_field(500, 0, "abc"), invalid_argument);                           // bad field index
}

BOOST_AUTO_TEST_CASE( set_individual_field_by_sample_type_mismatch ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  auto sample_names = vector<string>{ "NA12878", "NA12891", "NA12892" };
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  BOOST_CHECK_THROW(builder.set_integer_individual_field("AF", sample_names[0], {1, 2}), invalid_argument);
  BOOST_CHECK_THROW(builder.set_float_individual_field("ZIFMT", sample_names[0], {1.5f, 2.5f}), invalid_argument);
  BOOST_CHECK_THROW(builder.set_string_individual_field("AF", sample_names[0], "abc"), invalid_argument);
}

BOOST_AUTO_TEST_CASE( set_individual_fields_both_in_bulk_and_by_sample ) {
  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  auto bulk_changes = vector<int32_t>{1, 2, 3};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  builder.set_integer_individual_field("ZIFMT", bulk_changes);
  builder.set_integer_individual_field("ZIFMT", "NA12878", 5);

  // Attempting to build when both bulk changes and per-sample changes are queued up for a field should throw:
  BOOST_CHECK_THROW(builder.build(), logic_error);

  builder.clear().set_chromosome(0).set_alignment_start(5).set_ref_allele("A");
  builder.set_integer_individual_field("ZIFMT", bulk_changes);
  builder.set_integer_individual_field("GQ", "NA12878", 1);

  // Setting one field in bulk and another by sample should NOT throw, however:
  builder.build();
}


/*************************************************
 * Miscellaneous
 *************************************************/

BOOST_AUTO_TEST_CASE( disable_builder_validation ) {
  // WARNING: It's unsafe to test the "disable validation" feature in general, since attempting to do many
  //          of the actions that validation prevents will introduce undefined behavior into the test suite.
  //          However, it should be safe to at least do the check below with an invalid alignment stop.

  auto header = SingleVariantReader{"testdata/test_variants_for_variantbuilder.vcf"}.header();
  auto builder = VariantBuilder{header};
  builder.set_chromosome(0).set_alignment_start(5).set_ref_allele("A");

  // Disabling validation should turn off checks such as alignment stop < alignment start
  BOOST_CHECK_THROW(builder.set_alignment_stop(4).build(), logic_error);
  builder.set_enable_validation(false);
  builder.set_alignment_stop(4).build();
}
