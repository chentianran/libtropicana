#ifndef _TABLE_HH_
#define _TABLE_HH_

#include <unordered_map>
#include "key.hh"

class State;

struct KeyTable : public std::unordered_map<Key,State*>
{
    typedef std::unordered_map<Key,State*>::iterator Iterator;

    struct Insertion : public std::pair<Iterator,bool>
    {
        bool    inserted () const { return second; }
        bool    exists   () const { return ! second; }
        Key     key      () const { return first->first; }
        State*& state    () const { return first->second; }

        Insertion () { }
        Insertion (const std::pair<Iterator,bool>& i) : std::pair<Iterator,bool>(i) { }
    };

    Insertion try_insert (const Key& key)
    {
        return insert ( {key, 0} );
    }
};

#endif
