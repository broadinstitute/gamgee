#include "utils/short_value_optimized_storage.h"

#include <boost/test/unit_test.hpp>
#include <stdexcept>

using namespace std;
using namespace gamgee::utils;

void test_short_value_storage( const uint32_t capacity, const uint32_t short_value_upper_bound, const vector<vector<int32_t>>& values ) {
  auto storage = ShortValueOptimizedStorage<int32_t>{capacity, short_value_upper_bound};

  BOOST_CHECK_EQUAL(storage.capacity(), capacity);
  BOOST_CHECK_EQUAL(storage.num_values(), 0);
  BOOST_CHECK_EQUAL(storage.max_value_length(), 0);

  auto max_length = 0u;
  for ( auto i = 0u; i < storage.capacity(); ++i ) {
    BOOST_CHECK(! storage.is_set(i));
    storage.set(i, &(values[i][0]), values[i].size());
    BOOST_CHECK(storage.is_set(i));

    BOOST_CHECK_EQUAL(storage.capacity(), capacity);   // make sure capacity stays constant as we insert items
    BOOST_CHECK_EQUAL(storage.num_values(), i + 1);
    BOOST_CHECK_EQUAL(storage.value_length(i), values[i].size());

    if ( values[i].size() > max_length ) {
      max_length = values[i].size();
    }
    BOOST_CHECK_EQUAL(storage.max_value_length(), max_length);
  }

  for ( auto i = 0u; i < storage.capacity(); ++i ) {
    auto current_val = storage.get(i);
    BOOST_CHECK(current_val.first != nullptr);
    BOOST_CHECK_EQUAL(current_val.second, values[i].size());
    for ( auto j = 0u; j < current_val.second; ++j ) {
      BOOST_CHECK_EQUAL(current_val.first[j], values[i][j]);
    }
  }

  storage.clear();
  BOOST_CHECK_EQUAL(storage.capacity(), capacity);
  BOOST_CHECK_EQUAL(storage.num_values(), 0);
  BOOST_CHECK_EQUAL(storage.max_value_length(), 0);
}

BOOST_AUTO_TEST_CASE( test_all_short_values ) {
  auto short_values = vector<vector<int32_t>>{ {1}, {2, 3}, {4, 5, 6}, {7, 8} };

  test_short_value_storage(4, 5, short_values);  // capacity of 4 arrays, short value defined as <= length 5
}

BOOST_AUTO_TEST_CASE( test_all_long_values ) {
  auto long_values = vector<vector<int32_t>>{ {1, 2, 3, 4}, {5, 6, 7, 8, 9}, {10, 11, 12, 13} };

  test_short_value_storage(3, 2, long_values);   // capacity of 3 arrays, short value defined as <= length 2
}

BOOST_AUTO_TEST_CASE( test_mixed_short_and_long_values ) {
  auto mixed_values = vector<vector<int32_t>>{ {1, 2}, {3, 4, 5, 6, 7}, {8}, {9, 10}, {11, 12, 13, 14, 15, 16} };

  test_short_value_storage(5, 3, mixed_values);  // capacity of 5 arrays, short value defined as <= length 3
}

BOOST_AUTO_TEST_CASE( test_borderline_length_values ) {
  auto borderline_values = vector<vector<int32_t>>{ {1, 2, 3}, {4, 5, 6, 7}, {8, 9, 10} };

  test_short_value_storage(3, 3, borderline_values);  // capacity of 3 arrays, short value defined as <= length 3
}

BOOST_AUTO_TEST_CASE( test_set_by_vector ) {
  auto storage = ShortValueOptimizedStorage<int32_t>{5, 3};  // capacity of 5 arrays, short value defined as <= length 3
  auto mixed_values = vector<vector<int32_t>>{ {1, 2}, {3, 4, 5, 6, 7}, {8}, {9, 10}, {11, 12, 13, 14, 15, 16} };

  BOOST_CHECK_EQUAL(storage.capacity(), 5);
  BOOST_CHECK_EQUAL(storage.num_values(), 0);
  BOOST_CHECK_EQUAL(storage.max_value_length(), 0);

  for ( auto i = 0u; i < storage.capacity(); ++i ) {
    BOOST_CHECK(! storage.is_set(i));
    storage.set(i, mixed_values[i]);  // set by vector instead of by pointer
    BOOST_CHECK(storage.is_set(i));
    BOOST_CHECK_EQUAL(storage.capacity(), 5);   // make sure capacity stays constant as we insert items
    BOOST_CHECK_EQUAL(storage.num_values(), i + 1);
    BOOST_CHECK_EQUAL(storage.value_length(i), mixed_values[i].size());
    BOOST_CHECK_EQUAL(storage.max_value_length(), i < 1 ? 2 : (i < 4 ? 5 : 6));
  }

  for ( auto i = 0u; i < storage.capacity(); ++i ) {
    auto current_val = storage.get(i);
    BOOST_CHECK(current_val.first != nullptr);
    BOOST_CHECK_EQUAL(current_val.second, mixed_values[i].size());
    for ( auto j = 0u; j < current_val.second; ++j ) {
      BOOST_CHECK_EQUAL(current_val.first[j], mixed_values[i][j]);
    }
  }

  storage.clear();
  BOOST_CHECK_EQUAL(storage.capacity(), 5);
  BOOST_CHECK_EQUAL(storage.num_values(), 0);
  BOOST_CHECK_EQUAL(storage.max_value_length(), 0);
}

BOOST_AUTO_TEST_CASE( test_set_empty_values ) {
  auto storage = ShortValueOptimizedStorage<int32_t>{4, 5};  // capacity of 4 arrays, short value defined as <= length 5
  auto values = vector<vector<int32_t>>{ {1}, {} };

  BOOST_CHECK_EQUAL(storage.capacity(), 4);
  BOOST_CHECK_EQUAL(storage.num_values(), 0);
  BOOST_CHECK_EQUAL(storage.max_value_length(), 0);

  // Attempting to set empty values should have no effect on the state of the storage
  storage.set(0, nullptr, 1);           // null pointer
  BOOST_CHECK(! storage.is_set(0));
  storage.set(1, &(values[0][0]), 0);   // non-null pointer with num items == 0
  BOOST_CHECK(! storage.is_set(1));
  storage.set(2, values[1]);            // empty vector
  BOOST_CHECK(! storage.is_set(2));

  BOOST_CHECK_EQUAL(storage.num_values(), 0);    // num values should be unchanged
  BOOST_CHECK_EQUAL(storage.max_value_length(), 0);

  for ( auto i = 0u; i < storage.capacity(); ++i ) {
    auto val = storage.get(i);
    BOOST_CHECK(val.first == nullptr);
    BOOST_CHECK_EQUAL(val.second, 0);
  }
}

BOOST_AUTO_TEST_CASE( test_clear_individual_elements ) {
  auto storage = ShortValueOptimizedStorage<int32_t>{4, 3};

  for ( auto i = 0u; i < storage.capacity(); ++i ) {
    storage.set(i, vector<int32_t>{1});
  }

  storage.clear(1);
  storage.clear(3);

  BOOST_CHECK(storage.is_set(0));
  BOOST_CHECK_EQUAL(storage.get(0).first[0], 1);

  BOOST_CHECK(! storage.is_set(1));

  BOOST_CHECK(storage.is_set(2));
  BOOST_CHECK_EQUAL(storage.get(2).first[0], 1);

  BOOST_CHECK(! storage.is_set(3));
}

BOOST_AUTO_TEST_CASE( test_invalid_short_value_upper_bound ) {
  // short value upper bound must be > 0
  BOOST_CHECK_THROW(auto storage = ShortValueOptimizedStorage<int32_t>(4, 0), invalid_argument);
}

BOOST_AUTO_TEST_CASE( test_bounds_checking ) {
  auto storage = ShortValueOptimizedStorage<int32_t>{3, 2};
  auto val = vector<int32_t>{1, 2, 3};

  BOOST_CHECK_THROW(storage.set(3, &(val[0]), 3), out_of_range);
  BOOST_CHECK_THROW(storage.set(3, val), out_of_range);

  BOOST_CHECK_THROW(storage.set(4, &(val[0]), 3), out_of_range);
  BOOST_CHECK_THROW(storage.set(4, val), out_of_range);

  BOOST_CHECK_THROW(storage.get(3), out_of_range);

  BOOST_CHECK_THROW(storage.value_length(3), out_of_range);
}
