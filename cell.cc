#include "cell.hh"

Cell::Cell (const State& state) : RowMatrix (state.n - 1, state.n - 1)
{
    LOG (1, "creating cell for " << state.tab.active());

    int n = state.n - 1;
    int k = 0;
    for (int i = 0; i < state.m; ++i) {
        if (state.tab.key[i]) {
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
    lu.compute (*this);
    return std::round(fabs(lu.determinant()));
}
