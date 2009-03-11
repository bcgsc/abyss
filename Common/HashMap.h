#ifndef HASHMAP_H
#define HASHMAP_H 1

#include "config.h"

#if HAVE_UNORDERED_MAP
# include <unordered_map>
# define hash_multimap std::unordered_multimap
#elif HAVE_TR1_UNORDERED_MAP
# include <tr1/unordered_map>
# define hash_multimap std::tr1::unordered_multimap
#elif HAVE_BACKWARD_HASH_MAP || HAVE_EXT_HASH_MAP || HAVE_HASH_MAP
# undef __DEPRECATED
# if HAVE_BACKWARDS_HASH_MAP
#  include <backward/hash_map>
# elif HAVE_EXT_HASH_MAP
#  include <ext/hash_map>
# elif HAVE_HASH_MAP
#  include <hash_map>
# endif
using __gnu_cxx::hash_multimap;
#else
#  error A hash map implementation is required.
#endif

#if HAVE_UNORDERED_MAP || HAVE_TR1_UNORDERED_MAP
#elif HAVE_BACKWARD_HASH_MAP || HAVE_EXT_HASH_MAP || HAVE_HASH_MAP
// Implement a string hash function so that a string can be used as
// a key in STL maps and set.
namespace __gnu_cxx {
	template<> struct hash<std::string> {
		size_t operator()(const std::string& s) const {
			return hash<const char*>()(s.c_str());
		}
	};
}
#endif

#endif
