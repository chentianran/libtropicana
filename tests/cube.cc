#include <boost/test/unit_test.hpp>
#include <trop.hh>

#include <iostream>

BOOST_AUTO_TEST_CASE(cube)
{
    const int L = 4;
    RowMatrix A ((L+1)*(L+1)*(L+1),3);
    int b = 0;
    for (int i = 0; i <= L; ++i) {
        for (int j = 0; j <= L; ++j) {
            for (int k = 0; k <= L; ++k) {
                A.row(b)(0) = i;
                A.row(b)(1) = j;
                A.row(b)(2) = k;
                ++ b;
            }
        }
    }

    Trop trop;
    trop.compute (A);
    BOOST_CHECK_EQUAL (trop.volume, 6*L*L*L);
}
