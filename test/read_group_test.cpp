#include <string>

#include <boost/test/unit_test.hpp>

#include "sam/read_group.h"


using namespace std;
using namespace gamgee;
using namespace boost;

BOOST_AUTO_TEST_CASE( rg_simple_record_test ) {
  auto record = string{"@RG\tID:12345\tPL:Illumina\tSM:#$%^"};

  auto rg = ReadGroup(record);
  BOOST_CHECK_EQUAL(rg.id, "12345");
  BOOST_CHECK_EQUAL(rg.sample, "#$%^");
  BOOST_CHECK_EQUAL(rg.platform, "Illumina");
  BOOST_CHECK_EQUAL(rg.center, "");

  std::stringstream rg_stream;
  rg_stream << rg;
  BOOST_CHECK_EQUAL(rg_stream.str(), record);
}

BOOST_AUTO_TEST_CASE( rg_strange_punctuation_test) {
  auto record = string{"@RG\tID:lalala:lalala\tPL:Spiffy new platform"};

  auto rg = ReadGroup(record);
  BOOST_CHECK_EQUAL(rg.id, "lalala:lalala");
  BOOST_CHECK_EQUAL(rg.platform, "Spiffy new platform");

  std::stringstream rg_stream;
  rg_stream << rg;
  BOOST_CHECK_EQUAL(rg_stream.str(), record);
}

BOOST_AUTO_TEST_CASE( rg_user_defined_field_test ) {
  auto record = string{"@RG\tID:1\tab:aaaa\txy:xxxx"};

  auto rg = ReadGroup(record);
  BOOST_CHECK_EQUAL(rg.extra_fields.size(), 2);
  BOOST_CHECK_EQUAL(rg.extra_fields["xy"], "xxxx");
  BOOST_CHECK_EQUAL(rg.extra_fields["ab"], "aaaa");

  std::stringstream rg_stream;
  rg_stream << rg;
  BOOST_CHECK_EQUAL(rg_stream.str(), record);
}
