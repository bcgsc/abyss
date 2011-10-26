#ifndef HASHMAP_H
#define HASHMAP_H 1

#include "config.h"

#if HAVE_UNORDERED_MAP
# include <unordered_map>
# define hash_map std::unordered_map
# define hash_multimap std::unordered_multimap
#elif HAVE_TR1_UNORDERED_MAP
# include <tr1/unordered_map>
# define hash_map std::tr1::unordered_map
# define hash_multimap std::tr1::unordered_multimap
#else
# include <boost/unordered_map.hpp>
# define hash_map boost::unordered_map
# define hash_multimap boost::unordered_multimap
#endif

#endif
