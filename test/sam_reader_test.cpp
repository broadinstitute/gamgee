#include <boost/test/unit_test.hpp>
#include "sam_reader.h"

using namespace std;
using namespace foghorn;

BOOST_AUTO_TEST_CASE( read_bam ) 
{
    SamReader<SamIterator> reader{"testdata/test_simple.bam"};  // todo -- add sam/cram/vcf/bcf tests as well
    int counter = 0;
    for (const auto& sam : reader) {
        cout << sam << endl;  // todo -- substitute this for a check on the data values as we implement the Sam api
        ++counter;
    }
    BOOST_CHECK_EQUAL(counter, 33);
}
