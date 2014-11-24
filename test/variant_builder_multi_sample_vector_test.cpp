#include "variant_builder_multi_sample_vector.h"
#include "variant.h"

#include <boost/test/unit_test.hpp>
#include <stdexcept>

using namespace std;
using namespace gamgee;

// Input and expected output for both test_integer_multi_sample_vector_multiple_samples and test_integer_multi_sample_vector_multiple_samples_bulk_setting
const static auto multi_sample_tests_input_values = vector<vector<int32_t>>{ {1, 2, 3}, {1, 2}, {1}, {} };
const static auto multi_sample_tests_expected_output = vector<int32_t>{1, 2, 3, 1, 2, bcf_int32_vector_end, 1, bcf_int32_vector_end, bcf_int32_vector_end, bcf_int32_missing, bcf_int32_vector_end, bcf_int32_vector_end};

BOOST_AUTO_TEST_CASE( test_integer_multi_sample_vector_multiple_samples ) {
  // 4 samples, max of 3 values per sample
  const auto num_samples = 4u;
  const auto max_values_per_sample = 3u;
  auto vec = VariantBuilderMultiSampleVector<int32_t>{num_samples, max_values_per_sample, bcf_int32_missing, bcf_int32_vector_end};

  for ( auto sample_index = 0u; sample_index < multi_sample_tests_input_values.size(); ++sample_index ) {
    // Set each value for each sample individually
    for ( auto value_index = 0u; value_index < multi_sample_tests_input_values[sample_index].size(); ++value_index ) {
      vec.set_sample_value(sample_index, value_index, multi_sample_tests_input_values[sample_index][value_index]);
    }
  }

  BOOST_CHECK_EQUAL(vec.num_samples(), num_samples);
  BOOST_CHECK_EQUAL(vec.max_values_per_sample(), max_values_per_sample);
  BOOST_CHECK_EQUAL(multi_sample_tests_expected_output.size(), vec.get_vector().size());
  auto output_vector = vec.get_vector();
  for ( auto i = 0u; i < multi_sample_tests_expected_output.size(); ++i ) {
    BOOST_CHECK_EQUAL(multi_sample_tests_expected_output[i], output_vector[i]);
  }
}

BOOST_AUTO_TEST_CASE( test_integer_multi_sample_vector_multiple_samples_bulk_setting ) {
  // 4 samples, max of 3 values per sample
  const auto num_samples = 4u;
  const auto max_values_per_sample = 3u;
  auto vec = VariantBuilderMultiSampleVector<int32_t>{num_samples, max_values_per_sample, bcf_int32_missing, bcf_int32_vector_end};

  for ( auto sample_index = 0u; sample_index < multi_sample_tests_input_values.size(); ++sample_index ) {
    // Set all values for each sample at once
    vec.set_sample_values(sample_index, multi_sample_tests_input_values[sample_index]);
  }

  BOOST_CHECK_EQUAL(vec.num_samples(), num_samples);
  BOOST_CHECK_EQUAL(vec.max_values_per_sample(), max_values_per_sample);
  BOOST_CHECK_EQUAL(multi_sample_tests_expected_output.size(), vec.get_vector().size());
  auto output_vector = vec.get_vector();
  for ( auto i = 0u; i < multi_sample_tests_expected_output.size(); ++i ) {
    BOOST_CHECK_EQUAL(multi_sample_tests_expected_output[i], output_vector[i]);
  }
}

BOOST_AUTO_TEST_CASE( test_integer_multi_sample_vector_one_sample ) {
  // 1 sample, max of 3 values per sample
  const auto num_samples = 1u;
  const auto max_values_per_sample = 3u;
  auto vec = VariantBuilderMultiSampleVector<int32_t>{num_samples, max_values_per_sample, bcf_int32_missing, bcf_int32_vector_end};
  auto input_values = vector<int32_t>{1, 2};
  auto expected_output = vector<int32_t>{1, 2, bcf_int32_vector_end};

  vec.set_sample_values(0, input_values);

  BOOST_CHECK_EQUAL(vec.num_samples(), num_samples);
  BOOST_CHECK_EQUAL(vec.max_values_per_sample(), max_values_per_sample);
  BOOST_CHECK_EQUAL(expected_output.size(), vec.get_vector().size());
  auto output_vector = vec.get_vector();
  for ( auto i = 0u; i < expected_output.size(); ++i ) {
    BOOST_CHECK_EQUAL(expected_output[i], output_vector[i]);
  }
}

BOOST_AUTO_TEST_CASE( test_integer_multi_sample_vector_single_valued_field ) {
  // 4 samples, max of 1 value per sample
  const auto num_samples = 4u;
  const auto max_values_per_sample = 1u;
  auto vec = VariantBuilderMultiSampleVector<int32_t>{num_samples, max_values_per_sample, bcf_int32_missing, bcf_int32_vector_end};
  auto input_values = vector<vector<int32_t>>{ {1}, {2}, {}, {3} };
  auto expected_output = vector<int32_t>{1, 2, bcf_int32_missing, 3};

  for ( auto sample_index = 0u; sample_index < input_values.size(); ++sample_index ) {
    // Set all values for each sample at once
    vec.set_sample_values(sample_index, input_values[sample_index]);
  }

  BOOST_CHECK_EQUAL(vec.num_samples(), num_samples);
  BOOST_CHECK_EQUAL(vec.max_values_per_sample(), max_values_per_sample);
  BOOST_CHECK_EQUAL(expected_output.size(), vec.get_vector().size());
  auto output_vector = vec.get_vector();
  for ( auto i = 0u; i < expected_output.size(); ++i ) {
    BOOST_CHECK_EQUAL(expected_output[i], output_vector[i]);
  }
}

BOOST_AUTO_TEST_CASE( test_float_multi_sample_vector_multiple_samples ) {
  auto float_missing = 0.0f; bcf_float_set_missing(float_missing);
  auto float_vector_end = 0.0f; bcf_float_set_vector_end(float_vector_end);

  // 4 samples, max of 3 values per sample
  const auto num_samples = 4u;
  const auto max_values_per_sample = 3u;
  auto vec = VariantBuilderMultiSampleVector<float>{num_samples, max_values_per_sample, float_missing, float_vector_end};
  auto input_values = vector<vector<float>>{ {1.5, 2.5, 3.5}, {1.5, 2.5}, {1.5}, {} };
  auto expected_output = vector<float>{1.5, 2.5, 3.5, 1.5, 2.5, float_vector_end, 1.5, float_vector_end, float_vector_end, float_missing, float_vector_end, float_vector_end};

  for ( auto sample_index = 0u; sample_index < input_values.size(); ++sample_index ) {
    // Set each value for each sample individually
    for ( auto value_index = 0u; value_index < input_values[sample_index].size(); ++value_index ) {
      vec.set_sample_value(sample_index, value_index, input_values[sample_index][value_index]);
    }
  }

  BOOST_CHECK_EQUAL(vec.num_samples(), num_samples);
  BOOST_CHECK_EQUAL(vec.max_values_per_sample(), max_values_per_sample);
  BOOST_CHECK_EQUAL(expected_output.size(), vec.get_vector().size());
  auto output_vector = vec.get_vector();
  for ( auto i = 0u; i < expected_output.size(); ++i ) {
    if ( bcf_float_is_missing(expected_output[i]) ) {
      BOOST_CHECK(bcf_float_is_missing(output_vector[i]));
    }
    else if ( bcf_float_is_vector_end(expected_output[i]) ) {
      BOOST_CHECK(bcf_float_is_vector_end(output_vector[i]));
    }
    else {
      BOOST_CHECK_EQUAL(expected_output[i], output_vector[i]);
    }
  }
}

BOOST_AUTO_TEST_CASE( test_multi_sample_vector_zero_samples_and_values ) {
  // Ensure that we don't crash when creating empty vectors with num_samples or max_values_per_sample equal to 0

  // 0 samples, max of 0 values per sample
  auto vec = VariantBuilderMultiSampleVector<int32_t>{0u, 0u, bcf_int32_missing, bcf_int32_vector_end};
  BOOST_CHECK_EQUAL(vec.num_samples(), 0u);
  BOOST_CHECK_EQUAL(vec.max_values_per_sample(), 0u);

  // 4 samples, max of 0 values per sample
  vec = VariantBuilderMultiSampleVector<int32_t>{4u, 0u, bcf_int32_missing, bcf_int32_vector_end};
  BOOST_CHECK_EQUAL(vec.num_samples(), 4u);
  BOOST_CHECK_EQUAL(vec.max_values_per_sample(), 0u);

  // 0 samples, max of 4 values per sample
  vec = VariantBuilderMultiSampleVector<int32_t>{0u, 4u, bcf_int32_missing, bcf_int32_vector_end};
  BOOST_CHECK_EQUAL(vec.num_samples(), 0u);
  BOOST_CHECK_EQUAL(vec.max_values_per_sample(), 4u);
}
