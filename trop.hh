#ifndef _TROP_HH_
#define _TROP_HH_

#include "basic.hh"

struct Trop
{
    VolType volume;

    void compute (const Eigen::Ref<RowMatrix>& supp, bool compute_vol = true);
};

#endif
