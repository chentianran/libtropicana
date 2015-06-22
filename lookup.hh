#ifndef _LOOKUP_HH_
#define _LOOKUP_HH_

#include "basic.hh"

typedef vector<bool> Key;

struct Lookup : public VectorI
{
    int dim;
    int inactive;
    int end;
    Key key;

    inline Lookup () : inactive(0), end(0) { }

    inline Lookup (int n, int m) : VectorI(n+m), dim(n), inactive(n), end(n+m)
    {
        for (int i = 0; i < n; ++i)
            (*this)(i) = -1;
        for (int i = 0; i < m; ++i)
            (*this)(n+i) = i;
        key.assign (m, false);
    }

    inline void activate (int out, int in)
    {
        assert (in >= inactive && in < end);
        assert (out >= 0 && out < inactive);
        assert (false == key[(*this)(in)]);
        assert ((*this)(out) < 0);

        key[ (*this)(in) ] = true;              // set the corresponding bit in the key to 1
        (*this)(out) = (*this)(in);             // move it into the active group
        (*this)(in)  = (*this)(end - 1);        // last inactive move into this slot
        -- end;                                 // we got one less inactive constraint
        #ifndef NDEBUG
            (*this)(end) = -2;
        #endif
    }

    inline void pivot (int out, int in)
    {
        assert (in >= inactive && in < end);    // the one coming in should be inactive
        assert (out >= 0 && out < inactive);    // the one going out should be active
        assert (false == key[(*this)(in)]);

        if ((*this)(out) < 0)
            activate (out, in);
        else {
            key[(*this)(in) ] = true;
            key[(*this)(out)] = false;
            std::swap ((*this)(in), (*this)(out));
        }
    }

    inline void remove (int key)
    {
        assert (key >= inactive && key < end);      // the one coming in should be inactive

        (*this)(key) = (*this)(end - 1);
        -- end;
    }

    inline VectorI::ConstSegmentReturnType active() const
    {
        return segment(0,inactive);
    }

    inline void pivot_from (const Lookup& tab, int out, int in)
    {
        (*this) = tab;
        pivot (out, in);
    }

    inline bool check() const
    {
        for (int i = 0; i < inactive; ++i)
            if ((*this)(i) >= 0)
                assert (true == key[(*this)(i)]);
        for (int i = inactive; i < end; ++i)
            assert (false == key[(*this)(i)]);
    }
};

#endif
