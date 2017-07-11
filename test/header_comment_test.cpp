#include <sstream>
#include <string>

#include <boost/test/unit_test.hpp>

#include "sam/header_comment.h"

using namespace std;
using namespace gamgee;
using namespace boost;

BOOST_AUTO_TEST_CASE( co_simple_records_test ) {
  auto record = string{"@CO\trandom comment."};

  auto co = SamHeaderComment(record);

  BOOST_CHECK_EQUAL(co.comment, "random comment.");

  std::stringstream co_stream;
  co_stream << co;
  BOOST_CHECK_EQUAL(co_stream.str(), record);
}

BOOST_AUTO_TEST_CASE( co_user_defined_field_test ) {
  auto record = string{"@CO\tcomment\twith :tab."};

  auto co = SamHeaderComment(record);
  BOOST_CHECK_EQUAL(co.comment, "comment\twith :tab.");

  std::stringstream co_stream;
  co_stream << co;
  BOOST_CHECK_EQUAL(co_stream.str(), record);
}

