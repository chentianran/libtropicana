#ifndef _INVERSE_HH_
#define _INVERSE_HH_

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
                col(j) = D.col(j) - (D.col(j).dot(v) * col(k));

    }
};

#endif
