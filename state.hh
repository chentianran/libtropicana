#ifndef _STATE_HH_
#define _STATE_HH_

#include <utility>
#include <Eigen/Core>
#include "inverse.hh"
#include "lookup.hh"
#include "core.hh"

/// A state object describes a feasible solution of an inequality Ax ?= b.
/// It also holds information useful in further computation.
class State
{
public:

    const Core&      A;         ///< The constraint matrix for the LP problem
    const ColVector& b;         ///< The right hand side for the LP problem
    int              n;         ///< Dimension of the LP problem (i.e. the n.o. columns in A and length of x)
    int              m;         ///< Total number of constraints (i.e. the n.o. rows in A)
    Inverse          inv;       ///< Basic inverse Matrix (inverse of the basic matrix)
    Lookup           tab;       ///< Lookup table for active and inactive contraints
    ColVector        x;         ///< The state variable x
    ColVector        res;       ///< The residual vector Ax - b;
    Key              known_dir; ///< The bit mask that keeps track of known direction (1 for known and 0 for unknown)

    /// The data structure that holds information useful for pivoting.
    /// It keeps track of the dual index (offset in lookup talbe as well as
    /// actual row index in A) of an incoming constraint together with the
    /// step size needed to reach this constraint.
    struct PivotInfo
    {
        int    tab_id;          ///< The offset in the lookup table of the incoming constraint.
        int    row_id;          ///< The actual row index (in A) of the incoming constraint
        double step;            ///< The step size
    };

    /// Constructor that initializes the system of inequalities.
    /// Note that a state only keeps the references to A and b and does not
    /// copy the data.
    State (const Core& _A, const ColVector& _b) :
        A(_A),
        b(_b),
        inv(_A.cols()),
        tab(A.cols(),A.rows()),
        x(_A.cols()),
        res(_A.rows())
    {
        assert (A.rows() == b.size());
        assert (A.rows() >= A.cols());
        n = _A.cols();
        m = _A.rows();
        known_dir.reset();

        #ifndef NDEBUG
            for (int i = 0; i < m; ++i)
                assert (A(i,n-1) == -1.0);
        #endif
    }

    /// Move a feasible solution toward a given direction.
    /// This function moves a feasible solution toward the k-th outgoing 
    /// direction given by the basic inverse matrix.
    ///
    /// \param  Ad  The temporary buffer for A*d. This can be used for
    ///             updating the residual later.
    /// \param  k   The index of the outgoing direction. That is we are
    ///             moving along the k-th column of the basic inverse matrix.
    /// \param  sgn The direction. -1 for the reversed direction.
    /// 
    /// \return     This function returns true if the move is possible or
    ///             false if the k-th direction points toward outside of the
    ///             feasible region or an unbounded direction.
    bool move (ColVector& Ad, int k, double sgn = 1.0);

    void phase1();

    void branch_out (ColMatrix& AD) const;

    /// First half of the computation for moving from one vertex to an adjacent one.
    /// Given an edge direction and the intermediate data produced by "branching",
    /// this function perform the first half of the computation for leaving a vertex.
    /// It determines the incoming constraints and the necessary step size to reach it.
    /// The result is kept in a PivotInfo structure.
    /// Note that this function does not alter the state itself.
    PivotInfo leave (const ColMatrix& AD, int k) const;

    /// Complete the process of moving from one vertex to an adjacent one started by "leave".
    /// This is the actual pivoting step in the whole process.
    /// Using the intermediate data created by "branching" and "leave",
    /// this function creates a new state that represents the new vertex
    /// reached by leaving the given vertex along the given edge.
    ///
    /// \param src    The source state.
    /// \param AD     The intermediate data A*D produced by "branching".
    /// \param k      The index of the edge direction along which we are leaving a vertex.
    /// \param info   The intermediate pivoting data produced by "try_leaving".
    void arrive (
            const State& src, 
            const ColMatrix& AD, 
            int k, 
            const PivotInfo& info);

    bool check() const;
};

#endif
