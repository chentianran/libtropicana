#ifndef _CELL_HH_
#define _CELL_HH_

#include "state.hh"

struct Cell : public RowMatrix
{
    VectorI indices;

    Cell (const State& state);
    VolType vol();
};

#endif
