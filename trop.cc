#define EIGEN_RUNTIME_NO_MALLOC

#include <unordered_set>
#include "trop.hh"
#include "state.hh"
#include "cell.hh"

void Trop::more_verbose() { ++ _trop_verbose; }
void Trop::less_verbose() { -- _trop_verbose; }

State* Trop::new_state (const RowMatrix& A, const ColVector& b)
{
    State* new_state;
    if (! heap.empty()) {
        new_state = heap.front();
        heap.pop_front();
    } else {
        Eigen::internal::set_is_malloc_allowed (true);
        new_state = new State (A,b);
        Eigen::internal::set_is_malloc_allowed (false);
    }

    return new_state;
}

void Trop::del_state (State* state)
{
    heap.push_front (state);
}

void Trop::compute (const Eigen::Ref<RowMatrix>& supp, bool compute_vol)
{
    const int n = supp.cols();
    const int m = supp.rows();
    const int N = n + 1;

    int dup_discover = 0;
    int max_pool = 0;
    int total_cells = 0;

    RowMatrix A  (m, n+1);              // support matrix
    ColMatrix AD (m, n+1);              // A * D buffer
    ColVector b  (m);                   // right hand side

    A.leftCols(n) = supp;
    A.col(n).fill (-1.0);
    b.setRandom();

    std::unordered_set<Key>   known;
    std::forward_list<State*> pool;

    State* state = new State(A,b);
    state->phase1();
    //state->conv();
    known.insert (state->tab.key);

    pool.push_front (state);

    volume = 0;

    long long pool_size = 1;

    Eigen::internal::set_is_malloc_allowed(false);

    while (! pool.empty()) {
        state = pool.front();
        pool.pop_front();

        state->branch_out (AD);

        for (int k = 0; k < N; ++k) {
            if (! state->known_dir[k]) {
                int out_row = state->tab(k);
                if (out_row >= 0) {
                    State::PivotInfo info = state->try_leave (AD, k);
                    if (info.tab_id >= 0) {
                        Key key = state->tab.key;
                        key.set   (info.row_id);
                        key.reset (out_row);
                        if (known.end() == known.find (key)) {
                            known.insert (key);
                            State* branch = new_state(A,b);
                            branch->finish_leave (*state, AD, k, info);
                            pool.push_front (branch);
                            ++ pool_size;

                            //cout << state->tab.active() << " -> " << branch->tab.active() << "  (leaving " << k << ")" << endl;

                        } else {
                            //cout << state->tab.active() << " -> " << info.row_id << " in  (leaving " << k << ")" << endl;
                            ++ dup_discover;
                        }
                    }
                }
            }
        }

        Cell cell (*state);

        volume += cell.vol();

        del_state (state);
        -- pool_size;

        ++ total_cells;

        LOG(1, "pool: " << pool_size);
    }

    LOG(0, "total cells: " << total_cells);
    LOG(0, "duplicated:  " << dup_discover);
    LOG(0, "trying:      " << _stat_try_leave);
    LOG(0, "finishing:   " << _stat_try_leave);
    LOG(0, "update inv:  " << _stat_update_inv);

    for (State* s : heap)
        delete s;
}
