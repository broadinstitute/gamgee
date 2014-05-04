#include "../gamgee/utils/utils.h"

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
