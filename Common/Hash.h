#ifndef HASH_H_
#define HASH_H_ 1

#include "config.h"

#if HAVE_STD_HASH
#ifdef __FUJITSU
# include <string>
#else
# include <functional>
#endif
using std::hash;
# define NAMESPACE_STD_HASH_BEGIN namespace std {
# define NAMESPACE_STD_HASH_END }
#elif HAVE_STD_TR1_HASH
# include <tr1/functional>
using std::tr1::hash;
# define NAMESPACE_STD_HASH_BEGIN namespace std { namespace tr1 {
# define NAMESPACE_STD_HASH_END } }
#else
# include <boost/functional/hash.hpp>
using boost::hash;
# define NAMESPACE_STD_HASH_BEGIN namespace boost {
# define NAMESPACE_STD_HASH_END }
#endif

#endif
