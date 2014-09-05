#include "../gamgee/utils/utils.h"
#include "../gamgee/zip.h"

#include <boost/test/unit_test.hpp>

using namespace gamgee::utils;

BOOST_AUTO_TEST_CASE( sequence_utils_reverse_complement_test ) 
{
  const auto seq = std::string{"TTGATCTCCGAT"};
  const auto rev = std::string{"ATCGGAGATCAA"};
  BOOST_CHECK_EQUAL(reverse_complement(seq), rev);                     // check that reversing only reverses the sequence (correctly)
  BOOST_CHECK_EQUAL(reverse_complement(reverse_complement(seq)), seq); // re-reverse should go back to original
  BOOST_CHECK_EQUAL(reverse_complement(reverse_complement("AGGTCGAGT")), "AGGTCGAGT");
  BOOST_CHECK_EQUAL(reverse_complement(reverse_complement("taaacgttaaaatatccctag")), "taaacgttaaaatatccctag");
  BOOST_CHECK_EQUAL(reverse_complement(reverse_complement("taNACTTTTTTTNNtccctag")), "taNACTTTTTTTNNtccctag");
}

BOOST_AUTO_TEST_CASE( sequence_utils_complement_test ) 
{
  auto inplace_seq = std::string{"AGCT"}; 
  const auto copy_seq = std::string{"AGCT"}; 
  const auto result = std::string{"TCGA"}; 
  BOOST_CHECK_EQUAL(complement(inplace_seq), result);
  BOOST_CHECK_EQUAL(inplace_seq, result);          // after taking the complement in place (above), this should be true.
  BOOST_CHECK_EQUAL(complement(copy_seq), result);
  BOOST_CHECK_NE(copy_seq, result);                // after taking the complement on const (above), the value of copy_seq shouldn't be changed
  BOOST_CHECK_EQUAL(complement("TTGATCTCCGAT"), "AACTAGAGGCTA");
  BOOST_CHECK_EQUAL(complement("AAAAAAAAAAAA"), "TTTTTTTTTTTT");                     
  BOOST_CHECK_EQUAL(complement("G"), "C");                     
}

BOOST_AUTO_TEST_CASE( zip_iterators_test ) {
  const auto v_odd = std::vector<uint32_t>{1, 3, 5, 7};
  const auto v_even = std::vector<uint32_t>{2, 4, 6, 8};
  const auto v_alpha = std::vector<char>{'A', 'B', 'C', 'D'};

  // simple test with two vectors of the same type
  auto i = 0u;
  for (const auto tup : zip(v_odd, v_even)) {
    BOOST_CHECK_EQUAL(tup.get<0>(), v_odd[i]);
    BOOST_CHECK_EQUAL(tup.get<1>(), v_even[i]);
    ++i;
  }

  // testing different types
  auto j = 0u;
  for (const auto tup : zip(v_odd, v_even, v_alpha)) {
    BOOST_CHECK_EQUAL(tup.get<0>(), v_odd[j]);
    BOOST_CHECK_EQUAL(tup.get<1>(), v_even[j]);
    BOOST_CHECK_EQUAL(tup.get<2>(), v_alpha[j]);
    ++j;
  }

  // testing boost::tie
  auto k = 0u;
  for (const auto tup : zip(v_odd, v_even, v_alpha)) {
    auto odd = 0u;
    auto even = 0u;
    auto alpha = 'Z';
    boost::tie(odd, even, alpha) = tup;
    BOOST_CHECK_EQUAL(odd, v_odd[k]);
    BOOST_CHECK_EQUAL(even, v_even[k]);
    BOOST_CHECK_EQUAL(alpha, v_alpha[k]);
    ++k;
  }
}
