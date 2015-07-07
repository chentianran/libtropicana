#include "state.hh"

bool State::leave (ColVector& Ad, int k, double sgn)
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

void State::phase1()
{
    ColVector Ad (m);

    x.setZero();
    x(n-1) = - b.maxCoeff() - 0.01;
    inv.setIdentity();

    res = A * x - b;

    for (int k = 0; k < n; ++k)
        if (! leave (Ad, k, 1.0))
            leave (Ad, k, -1.0);

    assert (check());
}

void State::branch_out (ColMatrix& AD) const
{
    for (int i = tab.inactive; i < tab.end; ++i) {
        int row_id = tab(i);
        AD.row(row_id) = A.row(row_id) * inv;
        int n_neg = 0;
        for (int k = 0; k < n; ++k) {
            if (AD.row(row_id)(k) < -1e-8)
                ++ n_neg;
        }
        if (0 == n_neg) {
            WARNING ("Redundant: " << row_id);
        }
    }
}

State::PivotInfo State::try_leave (const ColMatrix& AD, int k) const
{
    assert (k >= 0 && k < n);

    double min_step   = 0.0;
    int    min_tab_id = -1;
    int    min_row_id = -1;

    for (int i = tab.inactive; i < tab.end; ++i) {
        int row_id = tab(i);
        assert (row_id >= 0 && row_id < m);
        double AiDk = AD.col(k)(row_id);
        if (AiDk < -1e-8) {
            double step = - res(row_id) / AiDk;
            if (-1 == min_tab_id || step < min_step) {
                min_step = step;
                min_tab_id = i;
                min_row_id = row_id;
            }
        }
    }

    return PivotInfo {min_tab_id, min_row_id, min_step};
}

void State::finish_leave (const State& S, const ColMatrix& AD, int k, const PivotInfo& info)
{
    assert (info.tab_id >= S.tab.inactive && info.tab_id < S.tab.end);
    assert (info.row_id >= 0 && info.row_id < S.m);
    assert (info.step > 0.0);

    x = S.x + info.step * S.inv.col(k);
    inv.pivot_from (S.inv, k, A.row(info.row_id), AD.col(k)(info.row_id));
    tab.pivot_from (S.tab, k, info.tab_id);
    update_res_from (S.res, info.step);

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
