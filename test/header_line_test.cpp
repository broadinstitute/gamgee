#include <sstream>
#include <string>

#include <boost/test/unit_test.hpp>

#include "sam/header_line.h"


using namespace std;
using namespace gamgee;
using namespace boost;

BOOST_AUTO_TEST_CASE( hd_simple_records_test ) {
  // w/ sorting order field.
  auto record = string{"@HD\tVN:1.4\tSO:coordinate"};
  auto hd = SamHeaderLine(record);
  BOOST_CHECK_EQUAL(hd.version, "1.4");
  BOOST_CHECK(hd.sorting_order == SamHeaderLine::SortingOrder::COORDINAMTE);
  BOOST_CHECK(hd.grouping == SamHeaderLine::Grouping::NOT_SET);
  std::stringstream hd_stream;
  hd_stream << hd;
  BOOST_CHECK_EQUAL(hd_stream.str(), record);

  // w/ grouping field.
  record = string{"@HD\tVN:1.4\tGO:reference"};
  hd = SamHeaderLine(record);
  BOOST_CHECK_EQUAL(hd.version, "1.4");
  BOOST_CHECK(hd.sorting_order == SamHeaderLine::SortingOrder::NOT_SET);
  BOOST_CHECK(hd.grouping == SamHeaderLine::Grouping::REFERENCE);
  hd_stream.str("");
  hd_stream << hd;
  BOOST_CHECK_EQUAL(hd_stream.str(), record);
}

BOOST_AUTO_TEST_CASE( hd_user_defined_field_test ) {
  auto record = string{"@HD\tVN:1.4\tab:aaaa\txy:xxxx"};
  std::stringstream record_stream;

  auto hd = SamHeaderLine(record);
  BOOST_CHECK_EQUAL(hd.version, "1.4");
  BOOST_CHECK_EQUAL(hd.extra_fields.size(), 2);
  BOOST_CHECK_EQUAL(hd.extra_fields["xy"], "xxxx");
  BOOST_CHECK_EQUAL(hd.extra_fields["ab"], "aaaa");

  std::stringstream hd_stream;
  hd_stream << hd;
  BOOST_CHECK_EQUAL(hd_stream.str(), record);
}
