#ifndef _TROP_HH_
#define _TROP_HH_

#include <forward_list>
#include "basic.hh"

class State;

struct Trop
{
    static void more_verbose();
    static void less_verbose();

    VolType volume;

    std::forward_list<State*> heap;

    State* new_state (const RowMatrix& A, const ColVector& b);
    void del_state (State* state);

    void compute (const Eigen::Ref<RowMatrix>& supp, bool compute_vol = true);
};

#endif
