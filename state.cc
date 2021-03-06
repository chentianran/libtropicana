#include "state.hh"
//#include "reltab.hh"

bool State::move (ColVector& Ad, int k, double sgn)
{
    assert (k >= 0 && k < n);

    double min_step = 0.0;
    int    min_key  = -1;
    int    min_row  = -1;

    for (int i = tab.inactive; i < tab.end; ++i) {
        int row_id = tab(i);
        assert (row_id >= 0 && row_id < m);
        if ( sgn * (Ad(row_id) = A.row(row_id).dot(inv.col(k))) < -1e-8) {
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
    res = A*x - b;

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
        if (! move (Ad, k, 1.0))
            move (Ad, k, -1.0);

    assert (check());
}

void State::branch_out (ColMatrix& AD) const
{
    // TODO: skip redundant constraints
    // TODO: use relation table to skip more constraints

    for (int i = tab.inactive; i < tab.end; ++i) {
        int row_id = tab(i);
        for (int k = 0; k < n; ++k) {
            if (tab.active()(k) >= 0) {
                if (! known_dir[k]) {
                    double AiDk;
                    AiDk = A.row(row_id).dot (inv.col(k));
                    AD (row_id,k) = AiDk;
                }
            }
        }
//        AD.row(row_id) = A.row(row_id) * inv;
//        int n_neg = 0;
//        int k_neg = -1;
//        for (int k = 0; k < n; ++k) {
//            if (AD.row(row_id)(k) < -1e-8) {
//                ++ n_neg;
//               k_neg = k;
//            }
//        }
//        if (0 == n_neg) {
//            //WARNING ("Redundant: " << row_id);
//        }
//        if (1 == n_neg) {
//            int dir_id = tab.active()(k_neg);   // the active constraint
//            reltab.mutex (dir_id, row_id);
//        }
    }
}

State::PivotInfo State::leave (const ColMatrix& AD, int k) const
{
    COUNT(_stat_try_leave);

    assert (k >= 0 && k < n);

    double min_step   = 0.0;
    int    min_tab_id = -1;
    int    min_row_id = -1;

    for (int i = tab.inactive; i < tab.end; ++i) {      // for each constraint in the inactive group
        // TODO: skip redundant constraints
        // TODO: use relation table to skip more constraints
        int row_id = tab(i);                            // get the actual row index
        assert (row_id >= 0 && row_id < m);             // make sure the row index is valid
        double AiDk = AD.col(k)(row_id);                // get A[i] * D[k]
        if (AiDk < -1e-8) {                             // get A[i] * D[k] < 0
            double step = - res(row_id) / AiDk;         // step = -(A[i] - b[i]) / (A[i] * D[k])
            if (-1 == min_tab_id || step < min_step) {  // see if this step is even smaller than those before
                min_step = step;                        // keeps track of the smallest step size
                min_tab_id = i;                         // keeps track of the offset in lookup table
                min_row_id = row_id;                    // ...as well as the actual row index in A
                assert (min_step > 1e-8);
            }
        }
    }

    return PivotInfo {min_tab_id, min_row_id, min_step};
}

void State::arrive (const State& S, const ColMatrix& AD, int k, const PivotInfo& info)
{
    COUNT(_stat_arrive);

    assert (info.tab_id >= S.tab.inactive && info.tab_id < S.tab.end);
    assert (info.row_id >= 0 && info.row_id < S.m);
    assert (info.step > 0.0);

    x = S.x + info.step * S.inv.col(k);
    inv.pivot_from (S.inv, k, A.row(info.row_id), AD.col(k)(info.row_id));
    tab.pivot_from (S.tab, k, info.tab_id);

    for (int i = 0; i < m; ++i) {
        if (S.tab.key[i])
            res(i) = 0.0;
        else
            res(i) = S.res(i) + info.step * AD.col(k)(i);
    }
    res(info.row_id) = 0.0;
    res(S.tab.active()(k)) = info.step;

    known_dir.reset();
    known_dir.set(k);

    assert (check());
}

bool State::check() const
{
    assert (tab.check());

    //--- check residual -------------------------------------------------------
    ColVector r = A*x - b;
    for (int i = 0; i < n; ++i) {       // for each active constraint
        int row_id = tab(i);
        if (row_id >= 0)
        assert (fabs(r(row_id)) < 1e-8);
    }
    for (int i = tab.inactive; i < tab.end; ++i) {
        int row_id = tab(i);
        assert (fabs(r(row_id)) > 1e-8);
        assert (fabs(r(row_id) - res(row_id) < 1e-12));
    }

    //--- check inverse --------------------------------------------------------
    for (int k = 0; k < n; ++k) {
        if (tab(k) >= 0) {
            assert (tab(k) < m);
            RowVector ek = A.row(tab(k)) * inv;
            for (int j = 0; j < n; ++j)
                if (fabs(ek(j) - ( (j==k) ? 1.0 : 0.0)) > 1e-13)
                    ERROR("A[" << tab(k) << "] * D[" << j << "] = " << ek(j));
        }
    }

    return true;
}
