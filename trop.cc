#include <forward_list>
#include <unordered_set>
#include "trop.hh"
#include "state.hh"
#include "cell.hh"

void Trop::compute (const Eigen::Ref<RowMatrix>& supp, bool compute_vol)
{
    int n = supp.cols();
    int m = supp.rows();
    int N = n + 1;

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

            if (new_state->leave_from (*state, k)) {
                if (known.end() == known.find (new_state->tab.key)) {
                    known.insert (new_state->tab.key);
                    pool.push_front (new_state);
                } else {
                    heap.push_front (new_state);
                }
            } else {
                heap.push_front (new_state);
            }
        }

        LOGVAR (0, state->tab.active());
        Cell cell (*state);

        volume += cell.vol();

        heap.push_front (state);
    }
}
