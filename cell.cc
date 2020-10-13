#include "cell.hh"
#include <Eigen/LU>

Cell::Cell (const State& state) : RowMatrix (state.n - 1, state.n - 1), indices(state.n)
{
    LOG (1, "creating cell for " << state.tab.active());

    int n = state.n - 1;
    int k = 0;
    for (int i = 0; i < state.m; ++i) {                 // for each constraint 
        if (state.tab.key[i]) {                         // if it is active
            indices[k] = i;                             // save the index
            if (k < n)
                row(k) = state.A.row(i).head(n);
            else {
                assert (n == k);
                for (int j = 0; j < n; ++j)
                    row(j) -= state.A.row(i).head(n);
            }
            ++k;
        }
    }
    assert (n + 1 == k);
}

VolType Cell::vol()
{
    Eigen::PartialPivLU<RowMatrix> lu;
    lu.compute (*this);
    return std::round(fabs(lu.determinant()));
}
