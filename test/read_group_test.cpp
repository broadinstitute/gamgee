#include <string>

#include <boost/test/unit_test.hpp>

#include "read_group.h"


using namespace std;
using namespace gamgee;
using namespace boost;

BOOST_AUTO_TEST_CASE( simple_record_test ) {
  auto record = string{"@RG\tID:12345\tPL:Illumina\tSM:#$%^"};
  auto rg = ReadGroup(record);
  BOOST_CHECK_EQUAL(rg.id, "12345");
  BOOST_CHECK_EQUAL(rg.sample, "#$%^");
  BOOST_CHECK_EQUAL(rg.platform, "Illumina");
  BOOST_CHECK_EQUAL(rg.center, "");
}

BOOST_AUTO_TEST_CASE( strange_punctuation_test) {
  auto record = string{"@RG\tID:lalala:lalala\tPL:Spiffy new platform"};
  auto rg = ReadGroup(record);
  BOOST_CHECK_EQUAL(rg.id, "lalala:lalala");
  BOOST_CHECK_EQUAL(rg.platform, "Spiffy new platform");
}

BOOST_AUTO_TEST_CASE(user_defined_field_test) {
  auto record = string{"@RG\txy:xxxx\tab:aaaa"};
  auto rg = ReadGroup(record);
}
