#include "variant/variant.h"
#include "variant/synced_variant_reader.h"
#include "variant/synced_variant_iterator.h"
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

const auto synced_variant_vcf_inputs = vector<string>{"testdata/var_idx/test_variants_csi.vcf.gz", "testdata/var_idx/test_variants_tabix.vcf.gz"};
const auto synced_variant_bcf_inputs = vector<string>{"testdata/var_idx/test_variants.bcf"};

const auto synced_variant_chrom_full = "1,20,22";
const auto synced_variant_bp_full = "1:10000000-10000000,20:10001000-10001000,20:10002000-10002000,20:10003000-10003000,22:10004000-10004000";
const auto synced_variant_chrom_partial = "1";
const auto synced_variant_bp_partial = "20:10001000-10001000";

// one (different) record each
BOOST_AUTO_TEST_CASE( synced_variant_reader_partial_test ) {
  const auto intervals1 = synced_variant_chrom_partial;
  for (const auto input_files : {synced_variant_vcf_inputs, synced_variant_bcf_inputs}) {
    auto truth_index = 0u;
    const auto reader1 = SyncedVariantReader<SyncedVariantIterator>{input_files, intervals1};
    for (const auto& vec : reader1) {
      for (const auto& record : vec) {
        BOOST_CHECK_EQUAL(record.ref(), "T");
        BOOST_CHECK_EQUAL(record.chromosome(), 0u);
        BOOST_CHECK_EQUAL(record.alignment_start(), 10000000u);
        BOOST_CHECK_EQUAL(record.alignment_stop(), 10000000u);
        BOOST_CHECK_EQUAL(record.n_alleles(), 2u);
        BOOST_CHECK_EQUAL(record.n_samples(), 3u);
        BOOST_CHECK_EQUAL(record.id(), "db2342");
      }
      ++truth_index;
    }
    BOOST_CHECK_EQUAL(truth_index, 1u);
  }

  const auto intervals2 = synced_variant_bp_partial;
  for (const auto input_files : {synced_variant_vcf_inputs, synced_variant_bcf_inputs}) {
    auto truth_index = 0u;
    const auto reader2 = SyncedVariantReader<SyncedVariantIterator>{input_files, intervals2};
    for (const auto& vec : reader2) {
      for (const auto& record : vec) {
        BOOST_CHECK_EQUAL(record.ref(), "GG");
        BOOST_CHECK_EQUAL(record.chromosome(), 1u);
        BOOST_CHECK_EQUAL(record.alignment_start(), 10001000u);
        BOOST_CHECK_EQUAL(record.alignment_stop(), 10001001u);
        BOOST_CHECK_EQUAL(record.n_alleles(), 2u);
        BOOST_CHECK_EQUAL(record.n_samples(), 3u);
        BOOST_CHECK_EQUAL(record.id(), "rs837472");
      }
      ++truth_index;
    }
    BOOST_CHECK_EQUAL(truth_index, 1u);
  }
}
/*
1	10000000	db2342	T	C	80	PASS	AF=0.5;AN=6	GT:GQ:PL:AF	0/1:25:10,0,100:3.1,2.2	0/0:12:0,10,1000:3.1,2.2
20	10001000	rs837472	GG	AA	8.4	PASS	AF=0.5;AN=6	GT:GQ:PL:AF	0/1:35:10,0,100:3.1,2.2	0/0:35:0,10,100:3.1,2.2
22	4000	.	GAT	G,GATAT	.	PASS	AF=0.5,0;AN=6	GT:GQ:PL:AF	1/2:35:10,0,100,2,4,8:3.1,2.2	0/0:35:0,10,100,2,4,8:3.1,2.2
22	5000	.	GAT	G,GATAT	.	PASS	AF=0.5,.;AN=6	GT:GQ:PL:AF	0/1:35:10,0,100,.,4,.:3.1,.	0/0:.:0,10,100,2,4,8:3.1,2.2

22	10004000	.	GAT	G,GATAT	.	MISSED	AN=6	GT:GQ:PL:AS	0/0:35:0,10,100,2,4,8:ABA	1/1:35:10,100,0,2,4,8:ABA
22	10005000	.	GAT	G,GATAT	.	MISSED	AN=6	GT:GQ:PL:AS	0/0:.:0,10,100,2,4,8:ABA	1/1:35:10,100,0,2,4,8:.
*/
const auto synced_variant_sparse_inputs = vector<string>{"testdata/synced_sparse_1.vcf.gz", "testdata/synced_sparse_2.vcf.gz"};

const auto synced_variant_sparse_truth_present = vector<vector<bool>>{ {true, true}, {true, true}, {true, false}, {true, false}, {false, true}, {false, true} };
const auto synced_variant_sparse_truth_chrom = vector<string>{ "1", "20", "22", "22", "22", "22", "22" };
const auto synced_variant_sparse_truth_start = vector<uint>{ 10000000, 10001000, 4000, 5000, 10004000, 10005000 };

BOOST_AUTO_TEST_CASE( synced_variant_reader_sparse_test ) {
  auto pos_truth_index = 0u;
  const auto reader1 = SyncedVariantReader<SyncedVariantIterator>{synced_variant_sparse_inputs, synced_variant_chrom_full};
  for (const auto& vec : reader1) {
    // The vector will always be of the same size as the number of inputs
    BOOST_CHECK_EQUAL(vec.size(), synced_variant_sparse_inputs.size());
    auto record_index = 0u;
    for (const auto& record : vec) {
      BOOST_CHECK_EQUAL(!missing(record), synced_variant_sparse_truth_present[pos_truth_index][record_index]);
      if (synced_variant_sparse_truth_present[pos_truth_index][record_index]) {
        BOOST_CHECK_EQUAL(record.chromosome_name(), synced_variant_sparse_truth_chrom[pos_truth_index]);
        BOOST_CHECK_EQUAL(record.alignment_start(), synced_variant_sparse_truth_start[pos_truth_index]);
      }
      ++record_index;
    }
    ++pos_truth_index;
  }
  BOOST_CHECK_EQUAL(pos_truth_index, 6u);
}

BOOST_AUTO_TEST_CASE( synced_variant_reader_move_test ) {
  for (const auto input_files : {synced_variant_vcf_inputs, synced_variant_bcf_inputs}) {
    auto reader0 = SyncedVariantReader<SyncedVariantIterator>{input_files, synced_variant_chrom_full};
    auto reader1 = SyncedVariantReader<SyncedVariantIterator>{input_files, synced_variant_chrom_full};
    auto moved = check_move_constructor(reader1);

    auto record0 = reader0.begin().operator*();
    auto moved_record = moved.begin().operator*();

    BOOST_CHECK_EQUAL(record0[0].alignment_start(), moved_record[0].alignment_start());
  }
}

BOOST_AUTO_TEST_CASE( synced_variant_iterator_move_test ) {
  for (const auto input_files : {synced_variant_vcf_inputs, synced_variant_bcf_inputs}) {
    auto reader0 = SyncedVariantReader<SyncedVariantIterator>{input_files, synced_variant_chrom_full};
    auto iter0 = reader0.begin();
    auto reader1 = SyncedVariantReader<SyncedVariantIterator>{input_files, synced_variant_chrom_full};
    auto iter1 = reader1.begin();
    auto moved = check_move_constructor(iter1);

    auto record0 = *iter0;
    auto moved_record = *moved;

    BOOST_CHECK_EQUAL(record0[0].alignment_start(), moved_record[0].alignment_start());
  }
}

BOOST_AUTO_TEST_CASE( synced_variant_reader_nonexistent_file ) {
  // Single non-existent file
  BOOST_CHECK_THROW(SyncedVariantReader<SyncedVariantIterator>(vector<string>{"foo/bar/nonexistent.vcf"}, ""), FileOpenException);

  // Multiple files, one non-existent
  BOOST_CHECK_THROW(SyncedVariantReader<SyncedVariantIterator>(vector<string>{"testdata/var_idx/test_variants_csi.vcf.gz", "foo/bar/nonexistent.vcf"}, ""), FileOpenException);
}
