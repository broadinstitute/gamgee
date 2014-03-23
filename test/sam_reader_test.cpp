#include <boost/test/unit_test.hpp>
#include "sam_reader.h"

using namespace std;
using namespace foghorn;

BOOST_AUTO_TEST_CASE( read_bam ) 
{
    SamReader<SamIterator> reader{"testdata/test_simple.bam"};  // todo -- add sam/cram/vcf/bcf tests as well
    auto counter = 0u;
    for (const auto& sam : reader) {
        BOOST_CHECK_EQUAL(sam.name().substr(0, 15), "30PPJAAXX090125");
        BOOST_CHECK_EQUAL(sam.chromosome(), 0);
        ++counter;
    }
    BOOST_CHECK_EQUAL(counter, 33u);
}
