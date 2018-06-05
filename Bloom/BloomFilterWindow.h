#ifndef BLOOMFILTERWINDOW_H
#define BLOOMFILTERWINDOW_H 1

#include "Bloom.h"
#include "BloomFilter.h"
#include "Common/HashFunction.h"
#include "Common/Kmer.h"
#include "Common/IOUtil.h"
#include <algorithm>
#include <vector>
#include <iostream>

/**
 * A bloom filter that represents a window
 * within a larger bloom filter.
 */
class BloomFilterWindow : public Konnector::BloomFilter
{
public:

	/** Constructor. */
	BloomFilterWindow() : Konnector::BloomFilter() { };

	/** Constructor.
	 *
	 * @param fullBloomSize size in bits of the containing bloom filter
	 * @param startBitPos index of first bit in the window
	 * @param endBitPos index of last bit in the window
	 */
	BloomFilterWindow(size_t fullBloomSize, size_t startBitPos,
			size_t endBitPos, size_t hashSeed=0) :
		Konnector::BloomFilter(endBitPos - startBitPos + 1, hashSeed),
		m_fullBloomSize(fullBloomSize),
		m_startBitPos(startBitPos),
		m_endBitPos(endBitPos)
	{
		assert(startBitPos < fullBloomSize);
		assert(endBitPos < fullBloomSize);
		assert(startBitPos <= endBitPos);
	}

	/**
	 * Get the full size (in bits) of the bloom filter that
	 * this window is a part of.
	 */
	size_t fullBloomSize()
	{
		return m_fullBloomSize;
	}

	/** Get the start bit position for the window. */
	size_t startBitPos()
	{
		return m_startBitPos;
	}

	/** Get the end bit position for the window. */
	size_t endBitPos()
	{
		return m_endBitPos;
	}

	/** Return the size of the bit array. */
	size_t size() const
	{
		return Konnector::BloomFilter::size();
	}

	/** Return the number of elements with count >= max_count. */
	size_t popcount() const
	{
		return Konnector::BloomFilter::popcount();
	}

	/** Return the estimated false positive rate */
	double FPR() const
	{
		return Konnector::BloomFilter::FPR();
	}

	/** Return whether the specified bit is set. */
	bool operator[](size_t i) const
	{
		if (i >= m_startBitPos && i <= m_endBitPos)
			return Konnector::BloomFilter::operator[](i - m_startBitPos);
		return false;
	}

	/** Return whether the object is present in this set. */
	bool operator[](const Bloom::key_type& key) const
	{
		return (*this)[Bloom::hash(key, m_hashSeed) % m_fullBloomSize];
	}

	/** Add the object with the specified index to this set. */
	void insert(size_t i)
	{
		if (i >= m_startBitPos && i <= m_endBitPos)
			Konnector::BloomFilter::insert(i - m_startBitPos);
	}

	/** Add the object to this set. */
	void insert(const Bloom::key_type& key)
	{
		insert(Bloom::hash(key, m_hashSeed) % m_fullBloomSize);
	}

	/** Operator for reading a bloom filter from a stream. */
	friend std::istream& operator>>(std::istream& in, BloomFilterWindow& o)
	{
		o.read(in, BITWISE_OVERWRITE);
		return in;
	}

	/** Operator for writing the bloom filter to a stream. */
	friend std::ostream& operator<<(std::ostream& out, const BloomFilterWindow& o)
	{
		o.write(out);
		return out;
	}

	/** Read a bloom filter window from a stream. */
	void read(std::istream& in, BitwiseOp readOp = BITWISE_OVERWRITE)
	{
		Bloom::FileHeader header = Bloom::readHeader(in);
		assert(in);

		m_fullBloomSize = header.fullBloomSize;
		m_startBitPos = header.startBitPos;
		m_endBitPos = header.endBitPos;

		if (m_hashSeed != header.hashSeed) {
			if (readOp == BITWISE_OVERWRITE) {
				m_hashSeed = header.hashSeed;
			} else {
				std::cerr << "error: can't union/intersect bloom filters with "
					<< "different hash seed values\n";
				exit(EXIT_FAILURE);
			}
		}

		size_t bits = header.endBitPos - header.startBitPos + 1;

		if (m_size != bits) {
			if (readOp == BITWISE_OVERWRITE) {
				Konnector::BloomFilter::resize(bits);
			} else {
				std::cerr << "error: can't union/intersect bloom filters with "
					<< "different sizes\n";
				exit(EXIT_FAILURE);
			}
		}

		readBits(in, m_array, bits, 0, readOp);

		assert(in);
	}

	/** Write a bloom filter window to a stream. */
	void write(std::ostream& out) const
	{
		Bloom::FileHeader header;
		header.fullBloomSize = m_fullBloomSize;
		header.startBitPos = m_startBitPos;
		header.endBitPos = m_endBitPos;
		header.hashSeed = m_hashSeed;

		Bloom::writeHeader(out, header);
		assert(out);

		size_t windowSize = m_endBitPos - m_startBitPos + 1;
		out.write(m_array, (windowSize + 7)/8);
		assert(out);
	}

private:

	size_t m_fullBloomSize;
	size_t m_startBitPos, m_endBitPos;
};

#endif
