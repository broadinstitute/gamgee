#include <sstream>
#include <string>

#include <boost/test/unit_test.hpp>

#include "sam/program.h"


using namespace std;
using namespace gamgee;
using namespace boost;

BOOST_AUTO_TEST_CASE( pg_simple_record_test ) {
  auto record = string{"@PG\tID:12345\tCL:bwa-mem\tVN:v-1.0"};

  auto pg = Program(record);
  BOOST_CHECK_EQUAL(pg.id, "12345");
  BOOST_CHECK_EQUAL(pg.name, "");
  BOOST_CHECK_EQUAL(pg.version, "v-1.0");
  BOOST_CHECK_EQUAL(pg.command_line, "bwa-mem");

  std::stringstream pg_stream;
  pg_stream << pg;
  BOOST_CHECK_EQUAL(pg_stream.str(), record);
}

BOOST_AUTO_TEST_CASE( pg_strange_punctuation_test) {
  auto record = string{"@PG\tID:lalala la:lala\tCL:bwa-mem -k 10"};

  auto pg = Program(record);
  BOOST_CHECK_EQUAL(pg.id, "lalala la:lala");
  BOOST_CHECK_EQUAL(pg.command_line, "bwa-mem -k 10");

  std::stringstream pg_stream;
  pg_stream << pg;
  BOOST_CHECK_EQUAL(pg_stream.str(), record);
}

BOOST_AUTO_TEST_CASE( pg_user_defined_field_test ) {
  auto record = string{"@PG\tID:1\tab:aaaa\txy:xxxx"};
  std::stringstream record_stream;

  auto pg = Program(record);
  BOOST_CHECK_EQUAL(pg.extra_fields.size(), 2);
  BOOST_CHECK_EQUAL(pg.extra_fields["xy"], "xxxx");
  BOOST_CHECK_EQUAL(pg.extra_fields["ab"], "aaaa");

  std::stringstream pg_stream;
  pg_stream << pg;
  BOOST_CHECK_EQUAL(pg_stream.str(), record);
}
