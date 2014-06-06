#ifndef HASHFUNCTION_H
#define HASHFUNCTION_H 1

#include "city.h"
#include <stddef.h>
#include <stdint.h>

static inline uint64_t hashmem(const void *p, size_t n)
{
	return CityHash64(static_cast<const char*>(p), n);
}

static inline uint64_t hashmem(const void *p, size_t n, size_t seed)
{
	return CityHash64WithSeed(static_cast<const char*>(p), n, seed);
}

#endif
