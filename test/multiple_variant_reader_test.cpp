#include "variant/variant.h"
#include "variant/variant_header_builder.h"
#include "variant/multiple_variant_reader.h"
#include "variant/multiple_variant_iterator.h"
#include "test_utils.h"

#include <boost/test/unit_test.hpp>
#include <unordered_set>

using namespace std;
using namespace gamgee;

// copied / migrated from variant_reader_test

BOOST_AUTO_TEST_CASE( multi_variant_reader_validation )
{
  const std::vector<std::string> filenames1{"testdata/test_variants.vcf", "testdata/extra_header.vcf"};
  const std::vector<std::string> filenames2{"testdata/test_variants.vcf", "testdata/missing_header.vcf"};

  // validate mismatched headers by default
  for (const auto filenames_v : {filenames1, filenames2})
    BOOST_CHECK_THROW(
      auto reader = MultipleVariantReader<MultipleVariantIterator>(filenames_v),
      HeaderCompatibilityException
    );

  // don't validate mismatched headers
  for (const auto filenames_v : {filenames1, filenames2})
    auto reader = MultipleVariantReader<MultipleVariantIterator>(filenames_v, false);
}

const auto multi_diff_truth_record_count      = vector<uint32_t>{4, 1, 1, 1, 2, 1, 1, 1};
const auto multi_diff_truth_chromosome        = vector<uint32_t>{0, 1, 1, 1, 1, 2, 2, 2};
const auto multi_diff_truth_alignment_starts  = vector<uint32_t>{10000000, 10001000, 10001999, 10002000, 10003000, 10004000, 10005000, 10006000};
const auto multi_diff_truth_ref               = vector<string>{"T", "GG", "TAGTGQA", "TAGTGQA", "A", "GAT", "GAT", "GAT"};
const auto multi_diff_truth_n_alleles         = vector<uint32_t>{2, 2, 2, 2, 2, 3, 3, 3};
const auto multi_diff_truth_id                = vector<string>{"db2342", "rs837472", ".", ".", ".", ".", ".", "."};
auto truth_file_indices_mult_alt = vector<unordered_multiset<uint32_t>> {{0,1,1,1},{0},{1},{0},{0,1},{0},{0},{0}};

BOOST_AUTO_TEST_CASE( multiple_variant_reader_difference_test ) {
  auto truth_index = 0u;
  const auto reader = MultipleVariantReader<MultipleVariantIterator>{{"testdata/test_variants.vcf", "testdata/test_variants_multiple_alt.vcf"}, false};
  for (const auto& vec : reader) {
    auto expected_file_indices = truth_file_indices_mult_alt[truth_index];
    BOOST_CHECK_EQUAL(vec.size(), multi_diff_truth_record_count[truth_index]);
    for (const auto& pair : vec) {
      const auto& record = pair.first;
      BOOST_CHECK_EQUAL(record.chromosome(), multi_diff_truth_chromosome[truth_index]);
      BOOST_CHECK_EQUAL(record.alignment_start(), multi_diff_truth_alignment_starts[truth_index]);
      BOOST_CHECK_EQUAL(record.ref(), multi_diff_truth_ref[truth_index]);
      BOOST_CHECK_EQUAL(record.n_alleles(), multi_diff_truth_n_alleles[truth_index]);
      BOOST_CHECK_EQUAL(record.n_samples(), 3u);
      BOOST_CHECK_EQUAL(record.id(), multi_diff_truth_id[truth_index]);

      auto find_result = expected_file_indices.find(pair.second);
      BOOST_CHECK(find_result != expected_file_indices.end());
      expected_file_indices.erase(find_result);
    }
    BOOST_CHECK(expected_file_indices.empty());   // check that we've seen and erased all expected
    ++truth_index;
  }
  BOOST_CHECK_EQUAL(truth_index, 8u);
}

void multiple_variant_reader_sample_test(const vector<string> samples, const bool include, const uint desired_samples) {
  auto filenames = vector<string>{"testdata/mvar_rdr/test_variants.vcf", "testdata/mvar_rdr/test_variants.bcf"};

  auto reader = MultipleVariantReader<MultipleVariantIterator>{filenames, false, samples, include};
  BOOST_CHECK_EQUAL(reader.combined_header().n_samples(), desired_samples);
}

BOOST_AUTO_TEST_CASE( multiple_variant_reader_sites_only )
{
  multiple_variant_reader_sample_test(vector<string>{}, true, 0); // exclude all samples (sites-only)
}

BOOST_AUTO_TEST_CASE( multiple_variant_reader_include_all_samples )
{
  multiple_variant_reader_sample_test(vector<string>{}, false, 6);  // include all samples by setting include == false and passing an empty list
}

BOOST_AUTO_TEST_CASE( multiple_variant_reader_including )
{
  multiple_variant_reader_sample_test(vector<string>{"NA12878"}, true, 1); // include only NA12878
  multiple_variant_reader_sample_test(vector<string>{"NA12878", "NA12892"}, true, 2); // include both these samples
}

BOOST_AUTO_TEST_CASE( multiple_variant_reader_excluding )
{
  multiple_variant_reader_sample_test(vector<string>{"NA12891"}, false, 5);  // exclude only NA12891
  multiple_variant_reader_sample_test(vector<string>{"NA12891", "NA12878_bcf"}, false, 4);  // exclude both these samples
}

auto truth_file_indices_mvr_hdr = vector<unordered_set<uint32_t>> {{0,1},{0,1},{0,1},{0,1},{0},{0},{1},{1}};

BOOST_AUTO_TEST_CASE( multiple_variant_reader_headers_test ) {
  const auto file1 = "testdata/mvr_hdr/test1.vcf";
  const auto file2 = "testdata/mvr_hdr/test2.vcf";

  auto header1 = SingleVariantReader{file1}.header();
  auto header2 = SingleVariantReader{file2}.header();
  auto combined_header = VariantHeaderBuilder{header1}.merge(header2).build();

  auto truth_index = 0u;
  const auto reader = MultipleVariantReader<MultipleVariantIterator>{{file1, file2}};
  BOOST_CHECK(reader.combined_header() != header1);
  BOOST_CHECK(reader.combined_header() != header2);
  BOOST_CHECK(reader.combined_header() == combined_header);
  for (const auto& vec : reader) {
    auto expected_file_indices = truth_file_indices_mvr_hdr[truth_index];
    for (const auto& pair : vec) {
      const auto& variant = pair.first;
      BOOST_CHECK(variant.header() != combined_header);
      // order is determined by priority queue - hard to predict
      BOOST_CHECK(variant.header() == header1 || variant.header() == header2);

      auto find_result = expected_file_indices.find(pair.second);
      BOOST_CHECK(find_result != expected_file_indices.end());
      expected_file_indices.erase(find_result);
    }
    ++truth_index;
    BOOST_CHECK(expected_file_indices.empty());   // check that we've seen and erased all expected
  }
}

const auto gvcf_truth_ref               = vector<string>{"T", "C", "GG"};
const auto gvcf_truth_chromosome        = vector<uint32_t>{0, 0, 1};
const auto gvcf_truth_alignment_starts  = vector<uint32_t>{10000000, 20000000, 10001000};
const auto gvcf_truth_alignment_stops   = vector<uint32_t>{10000000, 20000123, 10001001};
const auto gvcf_truth_n_alleles         = vector<uint32_t>{2, 2, 2};
const auto gvcf_truth_id                = vector<string>{"db2342", ".", "rs837472"};
auto truth_file_indices_gvcf = vector<unordered_set<uint32_t>> {{0,1},{0,1},{0,1}};

BOOST_AUTO_TEST_CASE( gvcf_test_multiple ) {
  auto truth_index = 0u;
  const auto reader = MultipleVariantReader<MultipleVariantIterator>{vector<string>{"testdata/test.g.vcf", "testdata/test.g.bcf"}};
  for (const auto& vec : reader) {
    auto expected_file_indices = truth_file_indices_gvcf[truth_index];
    for (const auto& pair : vec) {
      const auto& record = pair.first;
      BOOST_CHECK_EQUAL(record.ref(), gvcf_truth_ref[truth_index]);
      BOOST_CHECK_EQUAL(record.chromosome(), gvcf_truth_chromosome[truth_index]);
      BOOST_CHECK_EQUAL(record.alignment_start(), gvcf_truth_alignment_starts[truth_index]);
      BOOST_CHECK_EQUAL(record.alignment_stop(), gvcf_truth_alignment_stops[truth_index]);
      BOOST_CHECK_EQUAL(record.n_alleles(), gvcf_truth_n_alleles[truth_index]);
      BOOST_CHECK_EQUAL(record.n_samples(), 3u);
      BOOST_CHECK_EQUAL(record.id(), gvcf_truth_id[truth_index]);

      auto find_result = expected_file_indices.find(pair.second);
      BOOST_CHECK(find_result != expected_file_indices.end());
      expected_file_indices.erase(find_result);
    }
    ++truth_index;
    BOOST_CHECK(expected_file_indices.empty());   // check that we've seen and erased all expected
  }
}

BOOST_AUTO_TEST_CASE( multiple_variant_reader_nonexistent_file ) {
  // Single non-existent file
  BOOST_CHECK_THROW(MultipleVariantReader<MultipleVariantIterator>(vector<string>{"foo/bar/nonexistent.vcf"}), FileOpenException);

  // Multiple files, one non-existent
  BOOST_CHECK_THROW(MultipleVariantReader<MultipleVariantIterator>(vector<string>{"testdata/test_variants.vcf", "foo/bar/nonexistent.vcf"}), FileOpenException);
}
