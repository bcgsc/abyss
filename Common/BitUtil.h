#ifndef BITUTILS_H
#define BITUTILS_H 1

#include "config.h"
#include <cstdlib> // for exit
#include <iostream>
#include <stdint.h>

/** The return value of the CPUID instruction. */
struct CPUID { unsigned a, b, c, d; };

/** Return the result of the CPUID instruction
 * or -1 if it is not supported.
 */
static inline CPUID cpuid(unsigned op)
{
	CPUID x;
#if __GNUC__ && __x86_64__
	__asm__("cpuid" : "=a" (x.a), "=b" (x.b), "=c" (x.c), "=d" (x.d)
			: "a" (op));
	return x;
#else
	(void)op;
	x.a = x.b = x.c = x.d = static_cast<unsigned>(-1);
	return x;
#endif
}

/** Return whether this processor has the POPCNT instruction. */
static inline bool havePopcnt() { return cpuid(1).c & (1 << 23); }

/** Ensure that this machine has the popcount instruction. */
static inline void checkPopcnt()
{
#if ENABLE_POPCNT
	if (!havePopcnt()) {
		std::cerr <<
"error: " PACKAGE " has been compiled to use the POPCNT\n"
"instruction, which this machine lacks. Recompile using\n"
"configure --disable-popcnt to disable this feature.\n";
		exit(EXIT_FAILURE);
	}
#endif
}

/** Return the Hamming weight of x. */
static inline uint64_t popcount(uint64_t x)
{
#if ENABLE_POPCNT && __GNUC__ && __x86_64__
  __asm__("popcnt %1,%0" : "=r" (x) : "r" (x));
  return x;
#else
  x = (x & 0x5555555555555555ULL) +
    ((x >> 1) & 0x5555555555555555ULL);
  x = (x & 0x3333333333333333ULL) +
    ((x >> 2) & 0x3333333333333333ULL);
  x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
  x = x + (x >>  8);
  x = x + (x >> 16);
  x = x + (x >> 32);
  return x & 0x7FLLU;
#endif
}

#endif
