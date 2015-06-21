#include <Eigen/Core>

#include "basic.hh"
#include "debug.hh"

class Inverse : public ColMatrix
{
public:

    inline Inverse (int dim) : ColMatrix (dim, dim)
    {
        assert (dim > 0);
    }

    inline void pivot (int k, const RowVector& v)
    {
        pivot (k, v, v.dot(col(k)));
    }

    inline void pivot (int k, const RowVector& v, double hint)
    {
        assert (k >= 0 && k < cols());

        col(k) /= hint;
        for (int j = 0; j < cols(); ++j)
            if (j != k)
                col(j) -= col(j).dot(v) * col(k);
    }

    inline void pivot_from (const Inverse& D, int k, const RowVector& v, double hint)
    {
        assert (cols() == D.cols());
        assert (rows() == D.rows());
        assert (k >= 0 && k < cols());

        col(k) = D.col(k) / hint;
        for (int j = 0; j < cols(); ++j)
            if (j != k)
                col(j) = D.col(j) - (col(j).dot(v) * col(k));

    }
};

struct Lookup : public VectorI
{
    int dim;
    int inactive;
    int end;

    inline Lookup () : inactive(0), end(0) { }

    inline Lookup (int n, int m) : VectorI(n+m), dim(n), inactive(n), end(n+m)
    {
        for (int i = 0; i < n; ++i)
            (*this)(i) = -1;
        for (int i = 0; i < m; ++i)
            (*this)(n+i) = i;
    }

    inline void activate (int out, int in)
    {
        assert (in >= inactive && in < end);
        assert (out >= 0 && out < inactive);
        assert ((*this)(out) < 0);

        (*this)(out) = (*this)(in);
        (*this)(in)  = (*this)(end - 1);
        -- end;
        #ifndef NDEBUG
            (*this)(end) = -2;
        #endif
    }

    inline void pivot (int out, int in)
    {
        assert (in  >= inactive && in  < end);      // the one coming in should be inactive
        assert (out >= 0        && out < inactive); // the one going out should be active

        if ((*this)(out) < 0)
            activate (out, in);
        else
            std::swap ((*this)(in), (*this)(out));
    }

    inline void remove (int key)
    {
        assert (key >= inactive && key < end);      // the one coming in should be inactive

        (*this)(key) = (*this)(end - 1);
        -- end;
    }

    inline VectorI::ConstSegmentReturnType active() const
    {
        return segment(0,dim);
    }
};

class State
{
public:

    const RowMatrix& A;      // constraint matrix
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

    bool leave (int k, double sgn = 1.0)
    {
        assert (k >= 0 && k < n);

        double min_step = 0.0;
        int    min_key  = -1;
        int    min_row  = -1;

        for (int i = tab.inactive; i < tab.end; ++i) {
            int row_id = tab(i);
            assert (row_id >= 0 && row_id < m);
            if ( (Ad(row_id) = sgn * A.row(row_id) * inv.col(k)) < 0.0) {
                double step = - res(row_id) / Ad(row_id);
                if (-1 == min_key || step < min_step) {
                    min_step = step;
                    min_key  = i;
                    min_row  = row_id;
                }
            }
        }

        if (-1 == min_key)
            return false;

        x += min_step * sgn * inv.col(k);
        inv.pivot (k, A.row(min_row), Ad(min_row));
        tab.pivot (k, min_key);

        res = A * x - b;    // TODO update residual
        return true;
    }

    bool check() const
    {
        assert ((A*x-b - res).norm() < 1e-14);
        for (int k = 0; k < n; ++k) {
            if (tab(k) >= 0) {
                assert (tab(k) < m);
                RowVector ek = A.row(tab(k)) * inv;
                for (int j = 0; j < n; ++j)
                    assert (fabs(ek(j) - ( (j==k) ? 1.0 : 0.0)) < 1e-13);
            }
        }
    }

    void phase1()
    {
        x.setZero();
        x(n-1) = - b.minCoeff() - 0.01;
        inv.setIdentity();
        res = A * x - b;

        LOGVAR(0,res.transpose());

        LOGVAR(0,x.transpose());
        LOGVAR(0,inv);

        cout << "=== Phase I ===" << endl;

        for (int k = 0; k < n; ++k)
            if (! leave (k, 1.0))
                leave (k, -1.0);

        assert (check());
    }
};
