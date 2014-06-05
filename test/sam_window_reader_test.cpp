#include "sam_window_reader.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( single_window_readers ) {
  const auto& filename = "testdata/test_simple.bam";
  auto record_count = 0;
  std::vector<int32_t> actual_record_counts = {};

  /*
  Full contract TBD.
  > Do windows overlap-- or end at-- start and end of contig?
  > Assuming we're 1-based, does this code have any off-by-one errors?

  Currently expected windows in testdata/test_simple.bam with window 20k, step 10k:

  200
  255
  257
  10399
  10479

  10399
  10479
  22910
  22924
  27878
  27982
  29967

  22910
  22924
  27878
  27982
  29967
  30014
  34825
  34871

  30014
  34825
  34871
  44899
  44965
  46597
  46605

  44899
  44965
  46597
  46605
  52694
  52712
  57504
  57554
  59531
  59570

  52694
  52712
  57504
  57554
  59531
  59570
  60330
  60387

  60330
  60387

  89147
  89245
  89543
  89597

  89147
  89245
  89543
  89597
  94613
  94780
  97204
  97289

  94613
  94780
  97204
  97289
  */

  for (const auto& sam_window : SingleSamWindowReader{filename, 20000, 10000}) {
    int32_t window_record_count = 0;
    for (const auto& sam : sam_window) {
      window_record_count++;
    }
    actual_record_counts.push_back(window_record_count);
    record_count += window_record_count;
  }

  std::vector<int32_t> expected_record_counts = {5, 7, 8, 7, 10, 8, 2, 4, 8, 4};
  BOOST_CHECK_EQUAL_COLLECTIONS(
      actual_record_counts.begin(), actual_record_counts.end(),
      expected_record_counts.begin(), expected_record_counts.end());
  BOOST_CHECK_EQUAL(record_count, 63);
}
