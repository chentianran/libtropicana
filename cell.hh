#ifndef _CELL_HH_
#define _CELL_HH_

#include <Eigen/LU>
#include "state.hh"

struct Cell : public RowMatrix
{
    Eigen::PartialPivLU<RowMatrix> lu;

    Cell (const State& state);
    VolType vol();
};

#endif
