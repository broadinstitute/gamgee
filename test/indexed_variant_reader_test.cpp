#include "variant.h"
#include "indexed_variant_reader.h"
#include "indexed_variant_iterator.h"
#include "test_utils.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

// copied / migrated from variant_reader_test

/*
  NOTE: Updated test_variants.bcf and var_idx directory via-
    bcftools view testdata/test_variants.vcf -o testdata/var_idx/test_variants.bcf -O b
    bcftools view testdata/test_variants.vcf -o testdata/var_idx/test_variants_csi.vcf.gz -O z
    cp testdata/var_idx/test_variants_csi.vcf.gz testdata/var_idx/test_variants_tabix.vcf.gz
    bcftools index testdata/var_idx/test_variants.bcf
    bcftools index testdata/var_idx/test_variants_csi.vcf.gz
    bcftools index testdata/var_idx/test_variants_tabix.vcf.gz -t
    cp testdata/var_idx/test_variants.bcf testdata/test_variants.bcf
 */
const auto indexed_variant_vcf_inputs = vector<string>{"testdata/var_idx/test_variants_csi.vcf.gz", "testdata/var_idx/test_variants_tabix.vcf.gz"};
const auto indexed_variant_bcf_inputs = vector<string>{"testdata/var_idx/test_variants.bcf"};

const auto indexed_variant_chrom_full = vector<string> {"1", "20", "22"};
const auto indexed_variant_bp_full = vector<string> {"1:10000000-10000000", "20:10001000-10001000", "20:10002000-10002000", "20:10003000-10003000", "22:10004000-10004000"};
const auto indexed_variant_chrom_partial = vector<string> {"1"};
const auto indexed_variant_bp_partial = vector<string> {"20:10001000-10001000"};

// one (different) record each
BOOST_AUTO_TEST_CASE( indexed_variant_reader_partial_test ) {
  for (const auto filename : indexed_variant_bcf_inputs) {

    auto truth_index = 0u;
    const auto reader1 = IndexedVariantReader<IndexedVariantIterator>{filename, indexed_variant_chrom_partial};
    for (const auto& record : reader1) {
      BOOST_CHECK_EQUAL(record.ref(), "T");
      BOOST_CHECK_EQUAL(record.chromosome(), 0u);
      BOOST_CHECK_EQUAL(record.alignment_start(), 10000000u);
      BOOST_CHECK_EQUAL(record.alignment_stop(), 10000000u);
      BOOST_CHECK_EQUAL(record.n_alleles(), 2u);
      BOOST_CHECK_EQUAL(record.n_samples(), 3u);
      BOOST_CHECK_EQUAL(record.id(), "db2342");
      ++truth_index;
    }
    BOOST_CHECK_EQUAL(truth_index, 1u);

    truth_index = 0u;
    const auto reader2 = IndexedVariantReader<IndexedVariantIterator>{filename, indexed_variant_bp_partial};
    for (const auto& record : reader2) {
      BOOST_CHECK_EQUAL(record.ref(), "GG");
      BOOST_CHECK_EQUAL(record.chromosome(), 1u);
      BOOST_CHECK_EQUAL(record.alignment_start(), 10001000u);
      BOOST_CHECK_EQUAL(record.alignment_stop(), 10001001u);
      BOOST_CHECK_EQUAL(record.n_alleles(), 2u);
      BOOST_CHECK_EQUAL(record.n_samples(), 3u);
      BOOST_CHECK_EQUAL(record.id(), "rs837472");
      ++truth_index;
    }
    BOOST_CHECK_EQUAL(truth_index, 1u);
  }
}

BOOST_AUTO_TEST_CASE( indexed_variant_reader_move_test ) {
  for (const auto filename : indexed_variant_bcf_inputs) {
    auto reader0 = IndexedVariantReader<IndexedVariantIterator>{filename, indexed_variant_chrom_full};
    auto reader1 = IndexedVariantReader<IndexedVariantIterator>{filename, indexed_variant_chrom_full};
    auto moved = check_move_constructor(reader1);

    auto record0 = reader0.begin().operator*();
    auto moved_record = moved.begin().operator*();

    BOOST_CHECK_EQUAL(record0.alignment_start(), moved_record.alignment_start());
  }
}

BOOST_AUTO_TEST_CASE( indexed_variant_iterator_move_test ) {
  for (const auto filename : indexed_variant_bcf_inputs) {
    auto reader0 = IndexedVariantReader<IndexedVariantIterator>{filename, indexed_variant_chrom_full};
    auto iter0 = reader0.begin();
    auto reader1 = IndexedVariantReader<IndexedVariantIterator>{filename, indexed_variant_chrom_full};
    auto iter1 = reader1.begin();
    auto moved = check_move_constructor(iter1);

    auto record0 = *iter0;
    auto moved_record = *moved;

    BOOST_CHECK_EQUAL(record0.alignment_start(), moved_record.alignment_start());
  }
}


BOOST_AUTO_TEST_CASE( indexed_variant_reader_nonexistent_file ) {
  // VCF itself doesn't exist
  BOOST_CHECK_THROW(IndexedVariantReader<IndexedVariantIterator>("foo/bar/nonexistent.vcf", vector<string>{}), FileOpenException);

  // VCF exists, but no accompanying index file
  BOOST_CHECK_THROW(IndexedVariantReader<IndexedVariantIterator>("testdata/unindexed/test_unindexed.vcf", vector<string>{}), IndexLoadException);
}

