#ifndef HASH_H_
#define HASH_H_ 1

#include "config.h"

#if HAVE_TR1_FUNCTIONAL
# include <tr1/functional>
using std::tr1::hash;
#elif HAVE_FUNCTIONAL
# include <functional>
using std::hash;
#else
# include <boost/functional/hash.hpp>
using boost::hash;
#endif

#endif
