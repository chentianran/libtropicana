#include <boost/test/unit_test.hpp>
#include <trop.hh>

BOOST_AUTO_TEST_CASE(square)
{
    const int W = 2, H = 2;
    RowMatrix A ((W+1)*(H+1),2);
    for (int i = 0; i <= H; ++i) {
        for (int j = 0; j <= W; ++j) {
            A.row(i*(W+1)+j)(0) = i;
            A.row(i*(W+1)+j)(1) = j;
        }
    }

    Trop trop;
    trop.compute (A);
    BOOST_CHECK_EQUAL (trop.volume, W * H * 2);
}
