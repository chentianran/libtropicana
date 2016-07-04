#ifndef _CORE_HH_
#define _CORE_HH_

#include "basic.hh"

class CoreDense : public RowMatrix
{
public:

    CoreDense (const Eigen::Ref<RowMatrix>& points) : 
        RowMatrix (points.rows(), points.cols() + 1)
    {
        int n = points.cols();

        leftCols(n) = points;
        col(n).fill (-1.0);
    }
};

class CoreSparse : public SparseRowMatrix
{
public:

    CoreSparse (const Eigen::Ref<RowMatrix>& points) : 
        SparseRowMatrix (points.rows(), points.cols() + 1)
    {
        typedef Eigen::Triplet<double> T;
        std::vector<T> data;

        int m = points.rows();
        int n = points.cols();

        assert (m > 0);
        assert (n > 0);
        assert (m > n);

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j)
                if (fabs(points(i,j)) > 1e-8)
                    data.push_back (T(i,j,points(i,j)));
            data.push_back (T(i, n, -1.0));
        }

        setFromTriplets (data.begin(), data.end());
    }
};

#ifdef USE_SPARSE
typedef CoreSparse Core;
#else
typedef CoreDense Core;
#endif

#endif
