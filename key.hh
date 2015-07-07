#ifndef _KEY_HH_
#define _KEY_HH_

#ifndef KEY_BITS
#define KEY_BITS 256
#endif

#include <bitset>

typedef std::bitset<KEY_BITS> Key;

#endif
