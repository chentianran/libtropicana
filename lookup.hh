#ifndef _LOOKUP_HH_
#define _LOOKUP_HH_

#include "basic.hh"
#include "key.hh"

/// Constraints lookup table.
/// This structure keeps track of the list of active, inactive, and useless 
/// constraints both as lists and as a bit-mask.
/// Internally, it keeps an array that stores the indices of constraints 
/// in the form of | active | inactive | useless |. 
/// Activation of a constraint is done by moving the index from the inactive 
/// to the active area.
/// A special case is -1 may be listed in the "active" area to indicate the
/// use of Euclidean basis as artificial constraints.
/// A bit-mask also keeps the active constraints as one's.
/// The bit-mask and the lists are always kept in sync.
struct Lookup : public VectorI
{
    int dim;        ///< Dimension of the LP problem space
    int inactive;   ///< Start offset of the inactive indices
    int end;        ///< End offset of useful indices
    Key key;        ///< Bit-mask showing active indices as one

    inline Lookup () : inactive(0), end(0) { }

    /// Constructor that initialize an lookup table to the given dimension.
    /// The initial state is that all constraints are inactive, no constraint
    /// is useless, and n Euclidean directions are active artificial constraints.
    inline Lookup (int n, int m) : VectorI(n+m), dim(n), inactive(n), end(n+m)
    {
        for (int i = 0; i < n; ++i)
            (*this)(i) = -1;
        for (int i = 0; i < m; ++i)
            (*this)(n+i) = i;
        key.reset();
    }

    /// Activate a constraint.
    /// It replaces an Euclidean direction with a currently inactive constraint.
    /// Here, both arguments are offsets in the lookup table (not the actual
    /// indices of the constraints.
    /// 
    /// \param  out  The offset of the Euclidean direction to be replaced.
    /// \param  in   The offset in table that contains the index of the 
    ///              inactive constraint to come in.
    inline void activate (int out, int in)
    {
        assert (in >= inactive && in < end);    // "in"  should be in the inactive area
        assert (out >= 0 && out < inactive);    // "out" should be in the active area
        assert (false == key[(*this)(in)]);     // bit-mask should agree
        assert ((*this)(out) < 0);              // "out" direction should be an Euclidean direction

        key.set ((*this)(in));                  // set the corresponding bit in the key to 1
        (*this)(out) = (*this)(in);             // move it into the active group
        (*this)(in)  = (*this)(end - 1);        // last inactive move into this slot
        -- end;                                 // we got one less inactive constraint
        #ifndef NDEBUG                          // for debugging
            (*this)(end) = -2;                  // empty area of the lookup table should be
        #endif                                  // ...filled with -2
    }

    /// Swap a pair of active and inactive constraints.
    /// If the "out" direction is an Euclidean direction then it equivalent to
    /// just activating the "in" constract.
    /// 
    /// \param  out  The offset of the out direction.
    /// \param  in   The offset in table that contains the index of the 
    ///              inactive constraint to come in.
    inline void pivot (int out, int in)
    {
        assert (in >= inactive && in < end);    // the one coming in should be inactive
        assert (out >= 0 && out < inactive);    // the one going out should be active
        assert (false == key[(*this)(in)]);     // the one coming in should not be active

        if ((*this)(out) < 0)                   // if the out direction is an Euclidean direction
            activate (out, in);                 // then simply activate the constraint coming in
        else {                                  
            key.set   ((*this)(in));
            key.reset ((*this)(out));
            std::swap ((*this)(in), (*this)(out));
        }
    }

    /// Remove an inactive constraint.
    /// This function removes an inactive constraint from the lookup table.
    /// The designed use case is that if it is known that a constraint
    /// can never be active, then it can be safely removed from the lookup table.
    inline void remove (int id)
    {
        assert (id >= inactive && id < end);    // the one coming in should be inactive

        (*this)(id) = (*this)(end - 1);         // move the last inactive constraint here
        -- end;                                 // we have one less inactive constraint now
    }

    /// Access the indices of active constraints.
    /// This function returns a reference to a vector that contains all the active constraints.
    inline VectorI::ConstSegmentReturnType active() const
    {
        return segment(0,inactive);
    }

    /// Copy constructor followed by a pivot.
    inline void pivot_from (const Lookup& tab, int out, int in)
    {
        (*this) = tab;
        pivot (out, in);
    }

    /// Check the correctness of the lookup table.
    inline bool check() const
    {
        assert (inactive >= 0);
        assert (end >= inactive);

        for (int i = 0; i < inactive; ++i)          // for each index in the active group
            if ((*this)(i) >= 0)                    // if it is not an Euclidean direction
                assert (true == key[(*this)(i)]);   // make sure the bit-mask also marks it

        for (int i = inactive; i < end; ++i)        // for each index in the inactive group
            assert (false == key[(*this)(i)]);      // make sure the bit-mask has 0 for it

        return true;
    }
};

#endif
