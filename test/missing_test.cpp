#include <boost/test/unit_test.hpp>
#include <limits>
#include <cmath>

#include "missing.h"
#include "htslib/vcf.h"

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( detect_missing_float ) {
  float missing_value;
  bcf_float_set_missing(missing_value);
  BOOST_CHECK(missing(missing_value));
}

BOOST_AUTO_TEST_CASE( distinguish_missing_float_from_other_nan_values ) {
  // end of vector should not be confused with missing
  float end_of_vector_value;
  bcf_float_set_vector_end(end_of_vector_value);
  BOOST_CHECK(! missing(end_of_vector_value));

  // quiet NaN should not be confused with missing
  float quiet_nan = std::numeric_limits<float>::quiet_NaN();
  BOOST_CHECK(std::isnan(quiet_nan));
  BOOST_CHECK(! missing(quiet_nan));
}
