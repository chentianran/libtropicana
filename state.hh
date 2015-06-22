#ifndef _STATE_HH_
#define _STATE_HH_

#include <Eigen/Core>
#include "inverse.hh"
#include "lookup.hh"

class State
{
public:

    const RowMatrix& A;     // constraint matrix
    const ColVector& b;     // right hand side
    int              n;
    int              m;     // total number of constraints
    Inverse          inv;   // Basic inverse Matrix
    Lookup           tab;   // lookup table for contraints
    ColVector        x;     // state variable x
    ColVector        Ad;    // buffer for A*d
    ColVector        res;   // residual = Ax - b;

    State (const RowMatrix& _A, const ColVector& _b) :
        A(_A),
        b(_b),
        inv(_A.cols()),
        tab(A.cols(),A.rows()),
        x(_A.cols()),
        Ad(_A.rows()),
        res(_A.rows())
    {
        assert (A.rows() == b.size());
        assert (A.rows() >= A.cols());
        n = _A.cols();
        m = _A.rows();
        #ifndef NDEBUG
            for (int i = 0; i < m; ++i)
                assert (A(i,n-1) == -1.0);
        #endif
    }

    void update_res (double step)
    {
        res = A * x - b;
    }

    void update_res_from (const ColVector& r, double step)
    {
        res = A * x - b;
    }

    bool leave (int k, double sgn = 1.0);
    bool leave_from (const State& S, int k);

    void phase1();

    bool check() const;
};

#endif
