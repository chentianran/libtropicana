#include "state.hh"

bool State::leave (int k, double sgn)
{
    assert (k >= 0 && k < n);

    double min_step = 0.0;
    int    min_key  = -1;
    int    min_row  = -1;

    for (int i = tab.inactive; i < tab.end; ++i) {
        int row_id = tab(i);
        assert (row_id >= 0 && row_id < m);
        if ( sgn * (Ad(row_id) = A.row(row_id) * inv.col(k)) < -1e-8) {
            double step = - sgn* res(row_id) / Ad(row_id);
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
    update_res (min_step);

    assert (check());

    return true;
}

bool State::leave_from (const State& S, int k)
{
    assert (k >= 0 && k < n);

    double min_step = 0.0;
    int    min_key  = -1;
    int    min_row  = -1;
    int    n_neg    =  0;   // number of negative entries

    for (int i = S.tab.inactive; i < S.tab.end; ++i) {
        int row_id = S.tab(i);
        assert (row_id >= 0 && row_id < m);
        if ( (Ad(row_id) = A.row(row_id) * S.inv.col(k)) < -1e-8) {
            double step = - S.res(row_id) / Ad(row_id);
            if (-1 == min_key || step < min_step) {
                min_step = step;
                min_key  = i;
                min_row  = row_id;
            }
            ++ n_neg;
        }
    }

    if (1 == n_neg)
        WARNING ("interior point: " << min_row);

    if (-1 == min_key)
        return false;

    x = S.x + min_step * S.inv.col(k);
    inv.pivot_from (S.inv, k, A.row(min_row), Ad(min_row));
    tab.pivot_from (S.tab, k, min_key);
    update_res_from (S.res, min_step);

    assert (check());

    return true;
}

void State::phase1()
{
    x.setZero();
    x(n-1) = - b.maxCoeff() - 0.01;
    inv.setIdentity();

    res = A * x - b;

    for (int k = 0; k < n; ++k)
        if (! leave (k, 1.0))
            leave (k, -1.0);

    assert (check());
}

bool State::check() const
{
    assert (tab.check());
    assert ((A*x-b - res).norm() < 1e-14);
    for (int k = 0; k < n; ++k) {
        if (tab(k) >= 0) {
            assert (tab(k) < m);
            RowVector ek = A.row(tab(k)) * inv;
            for (int j = 0; j < n; ++j)
                if (fabs(ek(j) - ( (j==k) ? 1.0 : 0.0)) > 1e-13)
                    ERROR("A[" << tab(k) << "] * D[" << j << "] = " << ek(j));
            assert (fabs(res(tab(k))) < 1e-13);
        }
    }

    return true;
}
