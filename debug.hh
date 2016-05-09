#ifndef _DEBUG_HH_
#define _DEBUG_HH_

#include <iostream>
#include <cassert>
#include <iterator>
#include <algorithm>
#include <string>
#include <sstream>

using std::cout;
using std::clog;
using std::cerr;
using std::endl;
using std::ostream_iterator;
using std::copy;
using std::string;

extern int _trop_verbose;           // verbose level

#define ERROR(message) { \
    std::stringstream s; \
    s << message; \
    throw s.str().c_str(); \
}

#define WARNING(message) { clog << "WARNING: " << message << endl; }

#ifndef NDEBUG
#define LOGVAR(i,VAR) { if(i <= _trop_verbose) { clog << #VAR " = " << VAR << endl; } }
#define LOG(i,msg)    { if(i <= _trop_verbose) { clog << msg << endl; } }
#else
#define LOGVAR(i,VAR) {  }
#define LOG(i,msg)    {  }
#endif

extern long long _stat_try_leave;       // number of times starting a "leave"
extern long long _stat_finish_leave;    // number of times finishing a "leave"
extern long long _stat_update_inv;  // number of times updating inverse

#ifndef NDEBUG
    #define COUNT(x) { ++x; }
#else
    #define COUNT(x) { }
#endif

#endif
