#ifndef HASHFUNCTION_H
#define HASHFUNCTION_H 1

#include "city.h"
#include <stddef.h>
#include <stdint.h>

static inline uint64_t hashmem(const void *p, size_t n)
{
	return CityHash64(static_cast<const char*>(p), n);
}

#endif
