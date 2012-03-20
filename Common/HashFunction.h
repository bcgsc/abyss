#ifndef HASHFUNCTION_H
#define HASHFUNCTION_H 1

#include <stddef.h>
#include <stdint.h>

uint32_t hashlittle(const void *key, size_t length, uint32_t initval);

static inline uint32_t hashmem(const void *p, size_t n)
{
	return hashlittle(p, n, 0);
}

#endif
