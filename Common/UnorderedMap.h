#ifndef UNORDEREDMAP_H
#define UNORDEREDMAP_H 1

#include "config.h"
#include "Common/Hash.h"

#if HAVE_UNORDERED_MAP
# include <unordered_map>
using std::unordered_map;
using std::unordered_multimap;
#elif HAVE_TR1_UNORDERED_MAP
# include <tr1/unordered_map>
using std::tr1::unordered_map;
using std::tr1::unordered_multimap;
#else
# include <boost/unordered_map.hpp>
using boost::unordered_map;
using boost::unordered_multimap;
#endif

#endif
