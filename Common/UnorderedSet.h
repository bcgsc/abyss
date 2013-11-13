#ifndef UNORDEREDSET_H
#define UNORDEREDSET_H 1

#include "config.h"
#include "Common/Hash.h"

#if HAVE_UNORDERED_SET
# include <unordered_set>
using std::unordered_set;
using std::unordered_multiset;
#elif HAVE_TR1_UNORDERED_SET
# include <tr1/unordered_set>
using std::tr1::unordered_set;
using std::tr1::unordered_multiset;
#else
# include <boost/unordered_set.hpp>
using boost::unordered_set;
using boost::unordered_multiset;
#endif

#endif
