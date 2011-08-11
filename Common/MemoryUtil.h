#ifndef MEMORYUTIL_H
#define MEMORYUTIL_H 1

#include "config.h"
#include <fstream>
#include <unistd.h> // for getpagesize

/** Return the number of bytes used by the data and stack segments.
 * @return -1 on error
 */
static inline ssize_t getMemoryUsage()
{
#if HAVE_GETPAGESIZE
	std::ifstream in("/proc/self/statm");
	size_t size, resident, share, text, lib, data;
	return in >> size >> resident >> share >> text >> lib >> data
		? data * getpagesize() : -1;
#else
	return -1;
#endif
}

#endif
