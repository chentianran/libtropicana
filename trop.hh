#ifndef _TROP_HH_
#define _TROP_HH_

#include "basic.hh"

struct Trop
{
    static void more_verbose();
    static void less_verbose();

    VolType volume;

    void compute (const Eigen::Ref<RowMatrix>& supp, bool compute_vol = true);
};

#endif
