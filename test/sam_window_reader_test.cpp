#include "sam_window_reader.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( single_window_readers ) {
  const auto& filename = "testdata/test_simple.bam";
  auto window_count = 0;
  int32_t record_count = 0;

  /*
  Full contract TBD. Do windows overlap-- or end at-- start and end of contig?
  Currently expected windows in testdata/test_simple.bam with window 20k, step 10k:

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
    ++window_count;
    record_count += sam_window.m_size;
  }
  BOOST_CHECK_EQUAL(window_count, 10);
  BOOST_CHECK_EQUAL(record_count, 63);
}
