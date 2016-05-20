#define EIGEN_RUNTIME_NO_MALLOC

#include "trop.hh"
#include "state.hh"
#include "table.hh"
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

    RowMatrix A  (m, n+1);                              // support matrix
    ColMatrix AD (m, n+1);                              // A * D buffer
    ColVector b  (m);                                   // right hand side

    A.leftCols(n) = supp;
    A.col(n).fill (-1.0);
    b.setRandom();                                      // set random lifting values

    std::forward_list<State*> pool;                     // the pool of states to be explored
    KeyTable keytab;                                    // hash table containing all known states
    KeyTable::Insertion ins;

    State* v = new State(A,b);                          // create a root vertex
    v->phase1();                                        // solve the phase-one problem and reach a vertex
    //v->conv();
    pool.push_front (v);                                // this root vertex is to be explored
    ins = keytab.try_insert (v->tab.key);               // inser the key of the root vertex into key table
    assert (ins.inserted());                            // insertion must be successful
    ins.state() = v;

    volume = 0;

    long long pool_size = 1;


    Eigen::internal::set_is_malloc_allowed(false);

    while (! pool.empty()) {
        v = pool.front();
        pool.pop_front();

        v->branch_out (AD);                     // start to branch out

        for (int k = 0; k < N; ++k) {                   // for each possible edge direction to leave this vertex
            if (! v->known_dir[k]) {                    // if this edge direction is not pointing to a known vertex
                int out_row = v->tab(k);                // get the actual row index of the corresponding constraint
                if (out_row >= 0) {                     // if this direction is not an Euclidean direction
                    State::PivotInfo piv;               // for keeping track of intermediate data for pivoting
                    piv = v->leave (AD, k);     // try to leave this vertex along the k-th direction
                    if (piv.tab_id >= 0) {              // if there is an incoming constraint
                        Key key = v->tab.key;           // let's create the key of the new vertex
                        key.set   (piv.row_id);         // ...suppose the incoming row is in
                        key.reset (out_row);            // ...and the outgoing row is gone
                        ins = keytab.try_insert (key);  // try to insert the new key into the key-table
                        if (ins.inserted()) {           // if insertion is successful, ie., no such key existed
                            State* br = new_state(A,b); // then the branch state is unknown, we create it now
                            br->arrive (*v,AD,k,piv);   // reach the new vertex by finishing the pivoting process
                            ins.state() = br;           // save the state pointer into the key-table
                            pool.push_front (br);       // put the new vertex into the pool to be explored later
                            ++ pool_size;               // keep track of the pool size
                        } else {                        // if the adjacent vertex is already known
                            State* adj = ins.state();   // get the pointer to adjacent state
                            if (adj)                    // if that adjacent state is still to be explored
                                adj->known_dir[k] = 1;  // flag the direction as known
                            COUNT(dup_discover);        // keep track of the duplicated discovery
                        }
                    }
                }
            }
        }                                               // finished exploring this vertex

        //assert (keytab.has(v->tab.key));              // the current node must already be known
        keytab[v->tab.key] = 0;                         // reset the cached pointer because we are destroying v

        Cell cell (*v);                                 // this vertex will now be converted into a cell
        volume += cell.vol();                           // compute the volume and keep the total volume
        del_state (v);                                  // this vertex is now useless and can be released
        -- pool_size;                                   // keep track of the pool size
        ++ total_cells;                                 // keep track of the total cells

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
