#ifndef MEMUTILS_H
#define MEMUTILS_H

#include "Common/UnorderedMap.h"

template <class UnorderedMap>
static inline size_t
approxMemSize(const UnorderedMap& map)
{
	typedef typename UnorderedMap::value_type Entry;
	size_t filled_bucket_bytes = map.size() *
		(sizeof(Entry) + 3 * sizeof(Entry *));
	size_t empty_bucket_bytes = size_t((1.0 - map.load_factor()) *
		map.bucket_count() * sizeof(Entry *));
	return filled_bucket_bytes + empty_bucket_bytes;
}

#endif
