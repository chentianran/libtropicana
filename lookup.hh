#ifndef _LOOKUP_HH_
#define _LOOKUP_HH_

struct Lookup : public VectorI
{
    int dim;
    int inactive;
    int end;

    inline Lookup () : inactive(0), end(0) { }

    inline Lookup (int n, int m) : VectorI(n+m), dim(n), inactive(n), end(n+m)
    {
        for (int i = 0; i < n; ++i)
            (*this)(i) = -1;
        for (int i = 0; i < m; ++i)
            (*this)(n+i) = i;
    }

    inline void activate (int out, int in)
    {
        assert (in >= inactive && in < end);
        assert (out >= 0 && out < inactive);
        assert ((*this)(out) < 0);

        (*this)(out) = (*this)(in);
        (*this)(in)  = (*this)(end - 1);
        -- end;
        #ifndef NDEBUG
            (*this)(end) = -2;
        #endif
    }

    inline void pivot (int out, int in)
    {
        assert (in  >= inactive && in  < end);      // the one coming in should be inactive
        assert (out >= 0        && out < inactive); // the one going out should be active

        if ((*this)(out) < 0)
            activate (out, in);
        else
            std::swap ((*this)(in), (*this)(out));
    }

    inline void remove (int key)
    {
        assert (key >= inactive && key < end);      // the one coming in should be inactive

        (*this)(key) = (*this)(end - 1);
        -- end;
    }

    inline VectorI::ConstSegmentReturnType active() const
    {
        return segment(0,dim);
    }

    inline void pivot_from (const Lookup& tab, int out, int in)
    {
        (*this) = tab;
        pivot (out, in);
    }
};

#endif
