#include "basic.hh"
#include "state.hh"
#include "state_array.hh"

#include <iostream>
#include <boost/timer/timer.hpp>

using namespace std;

int main ()
{
    int n = 100;
    int iters = 1000;

    Inverse D (n);
    D.setIdentity();

    ArrayState Da (n);

    RowMatrix B (n,n);
    RowMatrix V (iters, n);

    B.setIdentity();
    V.setRandom();

    for (int i = 0; i < iters; ++i) {
        int r = i % n;
        B.row (r) = V.row(i);
    }

    {
        boost::timer::auto_cpu_timer timer ("Eigen: %w (%u)\n");
        for (int i = 0; i < iters; ++i) {
            int r = i % n;
            D.pivot (r, V.row(i));
        }
    }

    {
        boost::timer::auto_cpu_timer timer ("Array: %w (%u)\n");
        for (int i = 0; i < iters; ++i) {
            int r = i % n;
            Da.pivot (r, V.row(i).data());
        }
    }

    RowMatrix A (10,3);
    ColVector b (10);
    A.setRandom();
    A.col(2).fill(-1.0);
    b.setRandom();

    State state (A,b);

    cout << A << endl;
    cout << b << endl;

    state.phase1();

    LOGVAR (0, state.tab.active());

    state.leave(0);
    assert (state.check());

    LOGVAR (0, state.tab.active());

    state.leave(1);
    assert (state.check());

    LOGVAR (0, state.tab.active());

    //cout << "B*D = \n" << (B*D) << endl << endl;

    return 0;
}
