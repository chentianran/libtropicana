#include "state.hh"


bool State::leave_from (const State& S, int k)
{
    assert (k >= 0 && k < n);

    double min_step = 0.0;
    int    min_key  = -1;
    int    min_row  = -1;

    for (int i = S.tab.inactive; i < S.tab.end; ++i) {
        int row_id = S.tab(i);
        assert (row_id >= 0 && row_id < m);
        if ( (Ad(row_id) = A.row(row_id) * S.inv.col(k)) < 0.0) {
            double step = - S.res(row_id) / Ad(row_id);
            if (-1 == min_key || step < min_step) {
                min_step = step;
                min_key  = i;
                min_row  = row_id;
            }
        }
    }

    if (-1 == min_key)
        return false;

    x = S.x + min_step * S.inv.col(k);
    inv.pivot_from (S.inv, k, A.row(min_row), Ad(min_row));
    tab.pivot_from (S.tab, k, min_key);
    update_res_from (S.res, min_step);

    assert (check());
    
    return true;
}
