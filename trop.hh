#ifndef _TROP_HH_
#define _TROP_HH_

#include <forward_list>
#include "basic.hh"
#include "core.hh"
#include "cell.hh"

class State;

struct Trop
{
    bool compute_volume;
    bool keep_cells;

    Trop() : compute_volume(true), keep_cells(true) { }
    ~Trop();

    static void more_verbose();
    static void less_verbose();

    VolType volume;

    std::forward_list<Cell*>  cells;
    std::forward_list<State*> heap;

    State* new_state (const Core& A, const ColVector& b);
    void del_state (State* state);

    void compute (const Eigen::Ref<RowMatrix>& supp);
};

#endif
