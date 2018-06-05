/**
 * A Bloom filter
 * Copyright 2013 Shaun Jackman
 */
#ifndef BLOOMFILTER_H
#define BLOOMFILTER_H 1

#include "Bloom/Bloom.h"
#include "Common/Kmer.h"
#include "Common/IOUtil.h"
#include "Common/BitUtil.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <boost/dynamic_bitset.hpp>

/*
 * Put `BloomFilter` class in `Konnector` namespace to avoid collision with BTL
 * `BloomFilter` class of the same name.
 */
namespace Konnector {

/** A Bloom filter. */
class BloomFilter
{
  public:

	/** Constructor. */
	BloomFilter() : m_size(0), m_hashSeed(0), m_array(NULL) { }

	/** Constructor. */
	BloomFilter(size_t n, size_t hashSeed=0) : m_size(n),
		m_hashSeed(hashSeed)
	{
		m_array = new char[(n + 7)/8]();
	}

	~BloomFilter()
	{
		delete[] m_array;
	}

	/** Return the size of the bit array. */
	size_t size() const { return m_size; }

	/** Return the population count, i.e. the number of set bits. */
	size_t popcount() const
	{
		size_t count = 0;
		size_t bytes = (m_size + 7) / 8;
		size_t numInts = bytes / sizeof(uint64_t);
		size_t leftOverBytes = bytes % sizeof(uint64_t);
		uint64_t* intPtr = reinterpret_cast<uint64_t*>(m_array);
		for (size_t i = 0; i < numInts; i++) {
			count += ::popcount(intPtr[i]);
		}
		for (size_t i = (bytes - leftOverBytes)*8; i < m_size; i++) {
			if ((*this)[i])
				count++;
		}
		return count;
	}

	/** Return the estimated false positive rate */
	double FPR() const
	{
		return (double)popcount() / size();
	}

	/** Return whether the specified bit is set. */
	bool operator[](size_t i) const
	{
		assert(i < m_size);
		return m_array[i / 8] & 1 << (7 - i % 8);
	}

	/** Return whether the object is present in this set. */
	bool operator[](const Bloom::key_type& key) const
	{
		return (*this)[Bloom::hash(key, m_hashSeed) % m_size];
	}

	/** Add the object with the specified index to this set. */
	void insert(size_t i)
	{
		assert(i < m_size);
		m_array[i / 8] |= 1 << (7 - i % 8);
	}

	/** Add the object to this set. */
	void insert(const Bloom::key_type& key)
	{
		insert(Bloom::hash(key, m_hashSeed) % m_size);
	}

	/** Operator for reading a bloom filter from a stream. */
	friend std::istream& operator>>(std::istream& in, BloomFilter& o)
	{
		o.read(in, BITWISE_OVERWRITE);
		return in;
	}

	/** Operator for writing the bloom filter to a stream. */
	friend std::ostream& operator<<(std::ostream& out, const BloomFilter& o)
	{
		o.write(out);
		return out;
	}

	/** Read a bloom filter from a stream. */
	void read(std::istream& in, BitwiseOp readOp = BITWISE_OVERWRITE)
	{
		Bloom::FileHeader header = Bloom::readHeader(in);
		assert(in);

		if (m_hashSeed != header.hashSeed) {
			if (readOp == BITWISE_OVERWRITE) {
				m_hashSeed = header.hashSeed;
			} else {
				std::cerr << "error: can't union/intersect bloom filters with "
					<< "different hash seeds\n";
				exit(EXIT_FAILURE);
			}
		}

		if (m_size != header.fullBloomSize) {
			if (readOp == BITWISE_OVERWRITE) {
				resize(header.fullBloomSize);
			} else {
				std::cerr << "error: can't union/intersect bloom filters with "
					<< "different sizes\n";
				exit(EXIT_FAILURE);
			}
		}

		size_t bits = header.endBitPos - header.startBitPos + 1;
		readBits(in, m_array, bits, header.startBitPos, readOp);
		assert(in);
	}

	/** Write a bloom filter to a stream. */
	void write(std::ostream& out) const
	{
		Bloom::FileHeader header;
		header.fullBloomSize = m_size;
		header.startBitPos = 0;
		header.endBitPos = m_size - 1;
		header.hashSeed = m_hashSeed;

		Bloom::writeHeader(out, header);
		assert(out);

		out.write(m_array, (m_size + 7)/8);
		assert(out);
	}

	/** Resize the bloom filter (wipes the current data) */
	void resize(size_t size)
	{
		if (m_size > 0 && m_array != NULL)
			delete[] m_array;

		m_array = new char[(size + 7)/8]();
		m_size = size;
	}

  protected:

	size_t m_size;
	size_t m_hashSeed;
	char* m_array;
};

} // end Konnector namespace

#endif
