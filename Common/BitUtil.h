#ifndef BITUTILS_H
#define BITUTILS_H 1

#include "config.h"
#include <cstdlib> // for exit
#include <stdint.h>
#include <cassert>
#include <iostream>

enum BitwiseOp { BITWISE_OVERWRITE, BITWISE_OR, BITWISE_AND };

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

static const bool hasPopcnt = havePopcnt();

/** Return the Hamming weight of x. */
static inline uint64_t popcount(uint64_t x)
{
#if HAVE_POPCNT && __GNUC__ && __x86_64__
	if (hasPopcnt) {
		__asm__("popcnt %1,%0" : "=r" (x) : "r" (x));
		return x;
	}
#endif
	x = (x & 0x5555555555555555ULL) +
		((x >> 1) & 0x5555555555555555ULL);
	x = (x & 0x3333333333333333ULL) +
		((x >> 2) & 0x3333333333333333ULL);
	x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
	x = x + (x >>  8);
	x = x + (x >> 16);
	x = x + (x >> 32);
	return x & 0x7FLLU;
}

/**
 * Memory copy with bit-level resolution.
 *
 * @param src source data to copy (byte array)
 * @param dest destination of copied data (byte array)
 * @param bits number of bits to copy from src, starting from
 *        first bit of first byte of src
 * @param bitOffset bit offset into dest
 * @param op how to combine existing bits in dest with bits
 *        src (options: BITWISE_OVERWRITE, BITWISE_OR,
 *        BITWISE_AND)
 */
static inline void copyBits(char* src, char* dest, size_t bits,
	size_t bitOffset = 0, BitwiseOp op = BITWISE_OVERWRITE)
{
	size_t bytes = (bits + 7) / 8;
	size_t fullBytes = (bits % 8 == 0) ? bytes : bytes - 1;
	size_t byteOffset = bitOffset / 8;
	unsigned char shift = bitOffset % 8;
	unsigned char carryMask = 0xFF << (8 - shift);
	for (size_t i = 0; i < fullBytes; i++) {
		if (op == BITWISE_OVERWRITE) {
			dest[byteOffset + i] &= carryMask;
			dest[byteOffset + i + 1] &= ~carryMask;
		}
		if (op == BITWISE_AND) {
			dest[byteOffset + i] &= src[i] >> shift | carryMask;
			dest[byteOffset + i + 1] &= src[i] << (8 - shift) | ~carryMask;
		}
		else {
			dest[byteOffset + i] |= src[i] >> shift;
			dest[byteOffset + i + 1] |= src[i] << (8 - shift);
		}
	}
	if (fullBytes < bytes) {
		unsigned char bitsInLastByte = bits % 8;
		unsigned char lastByteMask = 0xFF << (8 - bitsInLastByte);
		unsigned char lastCarryMask = lastByteMask << (8 - shift);
		char lastByte = src[bytes - 1];
		lastByte &= lastByteMask;
		size_t lastByteIndex = byteOffset + bytes - 1;
		if (op == BITWISE_OVERWRITE)
			dest[lastByteIndex] &= ~(lastByteMask >> shift);
		if (op == BITWISE_AND)
			dest[lastByteIndex] &= lastByte >> shift | ~(lastByteMask >> shift);
		else
			dest[lastByteIndex] |= lastByte >> shift;
		if (lastCarryMask > 0) {
			if (op == BITWISE_OVERWRITE)
				dest[lastByteIndex + 1] &= ~lastCarryMask;
			if (op == BITWISE_AND)
				dest[lastByteIndex + 1] &= lastByte << (8 - shift) | ~lastCarryMask;
			else
				dest[lastByteIndex + 1] |= lastByte << (8 - shift);
		}
	}
}

/**
 * Read data from a stream with bit-level resolution.
 *
 * @param in input byte stream
 * @param dest destination of copied data (byte array)
 * @param bits number of bits to read from in, starting from
 *        first bit of first byte of in
 * @param bitOffset bit offset into dest
 */
static inline void readBits(std::istream& in, char* dest, size_t bits,
	size_t bitOffset = 0, BitwiseOp op = BITWISE_OVERWRITE)
{
	(void)readBits;

	const size_t IO_BUFFER_SIZE = 32 * 1024;
	size_t byteOffset = bitOffset / 8;
	unsigned char shift = bitOffset % 8;

	if (op == BITWISE_OVERWRITE && shift == 0) {
		// Simple byte-aligned copy.
		size_t bytes = (bits + 7) / 8;
		size_t fullBytes = (bits % 8 == 0) ? bytes : bytes - 1;
		in.read(dest + byteOffset, fullBytes);
		assert(in);
		if (fullBytes < bytes) {
			char lastByte;
			in.read(&lastByte, 1);
			assert(in);
			copyBits(&lastByte, dest + byteOffset + fullBytes, bits % 8);
		}
	} else {
		// Non-byte-aligned copy. A portion of each src byte is
		// carried over into the next dest byte.
		char buffer[IO_BUFFER_SIZE];
		for(size_t i = 0; i < bits;) {
			size_t bitsRead = std::min(IO_BUFFER_SIZE * 8, bits - i);
			size_t bytesRead = (bitsRead + 7)/8;
			in.read(buffer, bytesRead);
			assert(in);
			copyBits(buffer, dest + byteOffset, bitsRead, i + shift, op);
			i += bitsRead;
		}
	}
}

#endif
