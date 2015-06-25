#define EIGEN_RUNTIME_NO_MALLOC

#include <forward_list>
#include <unordered_set>
#include "trop.hh"
#include "state.hh"
#include "cell.hh"

void Trop::more_verbose() { ++ _trop_verbose; }
void Trop::less_verbose() { -- _trop_verbose; }

void Trop::compute (const Eigen::Ref<RowMatrix>& supp, bool compute_vol)
{
    const int n = supp.cols();
    const int m = supp.rows();
    const int N = n + 1;

    int dup_discover = 0;
    int max_pool = 0;
    int total_cells = 0;

    RowMatrix A (m, n+1);
    ColVector b (m);

    A.leftCols(n) = supp;
    A.col(n).fill (-1.0);
    b.setRandom();

    std::unordered_set<Key>   known;
    std::forward_list<State*> pool;
    std::forward_list<State*> heap;

    State* state = new State(A,b);
    state->phase1();
    known.insert (state->tab.key);

    pool.push_front (state);

    volume = 0;

    long long pool_size = 1;

    while (! pool.empty()) {
        state = pool.front();
        pool.pop_front();

        for (int k = 0; k < N; ++k) {

            State* new_state;
            if (! heap.empty()) {
                new_state = heap.front();
                heap.pop_front();
            } else {
                new_state = new State (A,b);
            }

            Eigen::internal::set_is_malloc_allowed(false);
            if (new_state->leave_from (*state, k)) {
                if (known.end() == known.find (new_state->tab.key)) {
                    known.insert (new_state->tab.key);
                    pool.push_front (new_state);
                    ++ pool_size;
                } else {
                    heap.push_front (new_state);
                    ++ dup_discover;
                }
            } else {
                heap.push_front (new_state);
            }
            Eigen::internal::set_is_malloc_allowed(true);
        }

        Cell cell (*state);

        volume += cell.vol();

        heap.push_front (state);
        -- pool_size;

        ++ total_cells;

        LOG(1, "pool: " << pool_size);
    }

    LOG(0, "total cells: " << total_cells);
    LOG(0, "duplicated:  " << dup_discover);
    
    for (State* s : heap)
        delete s;
}
