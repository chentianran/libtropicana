#include "basic.hh"
#include "state.hh"
#include "state_array.hh"

#include <iostream>
#include <boost/timer/timer.hpp>
#include <unordered_set>

using namespace std;

int main ()
{
    try {
        int n = 100;
        int iters = 1000;

        Inverse D (n);
        D.setIdentity();

        ArrayState Da (n);

        RowMatrix B (n,n);
        RowMatrix V (iters, n);

        B.setIdentity();
        V.setRandom();

        for (int i = 0; i < iters; ++i) {
            int r = i % n;
            B.row (r) = V.row(i);
        }

        {
            boost::timer::auto_cpu_timer timer ("Eigen: %w (%u)\n");
            for (int i = 0; i < iters; ++i) {
                int r = i % n;
                D.pivot (r, V.row(i));
            }
        }

        {
            boost::timer::auto_cpu_timer timer ("Array: %w (%u)\n");
            for (int i = 0; i < iters; ++i) {
                int r = i % n;
                Da.pivot (r, V.row(i).data());
            }
        }

        const int W = 3, H = 4;
        RowMatrix A (W * H,3);
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                A.row(i*W+j)(0) = i;
                A.row(i*W+j)(1) = j;
            }
        }
        ColVector b (W*H);
        A.col(2).fill(-1.0);
        b.setRandom();

        cout << A << endl;
        cout << b << endl;

        std::unordered_set<Key> known;
        std::vector<State*>     pool;

        State* state = new State(A,b);
        state->phase1();
        known.insert (state->tab.key);

        pool.push_back (state);

        while (! pool.empty()) {
            state = pool.back();
            pool.pop_back();

            for (int k = 0; k < 3; ++k) {
                State* new_state = new State (A,b);
                if (new_state->leave_from (*state, k)) {
                    if (known.end() == known.find (new_state->tab.key)) {
                        LOGVAR (0, new_state->tab.active());
                        known.insert (new_state->tab.key);
                        pool.push_back (new_state);
                    }
                }
            }

            delete state;
        }

        //cout << "B*D = \n" << (B*D) << endl << endl;
    } catch (const char* err) {
        cerr << err << endl;
        return 1;
    }
    return 0;
}
