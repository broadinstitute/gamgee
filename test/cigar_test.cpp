#include <boost/test/unit_test.hpp>

#include "sam/cigar.h"

using namespace std;
using namespace gamgee;


BOOST_AUTO_TEST_CASE( cigar_consumes_read_bases ) {
	BOOST_CHECK(Cigar::consumes_read_bases(CigarOperator::M));
	BOOST_CHECK(Cigar::consumes_read_bases(CigarOperator::EQ));
	BOOST_CHECK(Cigar::consumes_read_bases(CigarOperator::X));
	BOOST_CHECK(Cigar::consumes_read_bases(CigarOperator::I));
	BOOST_CHECK(!Cigar::consumes_read_bases(CigarOperator::D));
	BOOST_CHECK(Cigar::consumes_read_bases(CigarOperator::S));
	BOOST_CHECK(!Cigar::consumes_read_bases(CigarOperator::H));
	BOOST_CHECK(!Cigar::consumes_read_bases(CigarOperator::N));
	BOOST_CHECK(!Cigar::consumes_read_bases(CigarOperator::P));
	BOOST_CHECK(!Cigar::consumes_read_bases(CigarOperator::B));
}

BOOST_AUTO_TEST_CASE( cigar_consumes_reference_bases ) {
	BOOST_CHECK(Cigar::consumes_reference_bases(CigarOperator::M));
	BOOST_CHECK(Cigar::consumes_reference_bases(CigarOperator::EQ));
	BOOST_CHECK(Cigar::consumes_reference_bases(CigarOperator::X));
	BOOST_CHECK(!Cigar::consumes_reference_bases(CigarOperator::I));
	BOOST_CHECK(Cigar::consumes_reference_bases(CigarOperator::D));
	BOOST_CHECK(!Cigar::consumes_reference_bases(CigarOperator::S));
	BOOST_CHECK(!Cigar::consumes_reference_bases(CigarOperator::H));
	BOOST_CHECK(Cigar::consumes_reference_bases(CigarOperator::N));
	BOOST_CHECK(!Cigar::consumes_reference_bases(CigarOperator::P));
	BOOST_CHECK(!Cigar::consumes_reference_bases(CigarOperator::B));
}
