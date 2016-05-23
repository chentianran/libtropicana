#ifndef _CORE_HH_
#define _CORE_HH_

#include "basic.hh"

class CoreDense : public RowMatrix
{
public:

    using RowMatrix::RowMatrix;
};

typedef CoreDense Core;

#endif
