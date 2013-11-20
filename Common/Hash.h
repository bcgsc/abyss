#ifndef HASH_H_
#define HASH_H_ 1

#include "config.h"

#if HAVE_STD_HASH
# include <functional>
using std::hash;
# define NAMESPACE_STD_HASH_BEGIN namespace std {
# define NAMESPACE_STD_HASH_END }
#elif HAVE_STD_TR1_HASH
# include <tr1/functional>
using std::tr1::hash;
# define NAMESPACE_STD_HASH_BEGIN namespace std { namespace tr1 {
# define NAMESPACE_STD_HASH_END } }
#else
# error One of either std::hash or std::tr1::hash is required.
#endif

#endif
