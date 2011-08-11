#ifndef MEMORYUTIL_H
#define MEMORYUTIL_H 1

#include <fstream>
#include <unistd.h> // for getpagesize

/** Return the number of bytes used by the data and stack segments.
 * @return -1 on error
 */
static inline ssize_t getMemoryUsage()
{
	std::ifstream in("/proc/self/statm");
	size_t size, resident, share, text, lib, data;
	return in >> size >> resident >> share >> text >> lib >> data
		? data * getpagesize() : -1;
}

#endif
