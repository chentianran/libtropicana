#include <boost/test/unit_test.hpp>
#include <trop.hh>

//#include <iostream>

BOOST_AUTO_TEST_CASE(hypercube)
{
    for (int D = 3; D < 8; ++D) {
        int D2 = 1;     // 2^D
        int Df = 1;     // D!
        for (int i = 0; i < D; ++i) {
            D2 *= 2;
            Df *= (i+1);
        }

        RowMatrix A (D2, D);            
        for (int i = 0; i < D2; ++i)
            for (int j = 0; j < D; ++j)
                A.row(i)(j) = (i & (1 << j)) ? 1.0 : 0.0;

        Trop trop;
        trop.compute (A);
        BOOST_CHECK_EQUAL (trop.volume, Df);
    }
}
