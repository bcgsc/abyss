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
class BloomFilterWindow : public BloomFilter
{
public:

	/** Constructor. */
	BloomFilterWindow() : BloomFilter() { };

	/** Constructor.
	 *
	 * @param fullBloomSize size in bits of the containing bloom filter
	 * @param startBitPos index of first bit in the window
	 * @param endBitPos index of last bit in the window
	 */
	BloomFilterWindow(size_t fullBloomSize, size_t startBitPos, size_t endBitPos) :
		BloomFilter(endBitPos - startBitPos + 1),
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
		return BloomFilter::size();
	}

	/** Return the number of elements with count >= MAX_COUNT. */
	size_t popcount() const
	{
		return BloomFilter::popcount();
	}

	/** Return the estimated false positive rate */
	double FPR() const
	{
		return BloomFilter::FPR();
	}

	/** Return whether the specified bit is set. */
	bool operator[](size_t i) const
	{
		if (i >= m_startBitPos && i <= m_endBitPos)
			return BloomFilter::operator[](i - m_startBitPos);
		return false;
	}

	/** Return whether the object is present in this set. */
	bool operator[](const Bloom::key_type& key) const
	{
		return (*this)[Bloom::hash(key) % m_fullBloomSize];
	}

	/** Add the object with the specified index to this set. */
	void insert(size_t i)
	{
		set(i, true);
	}

	/**
	 * Set or clear position i of the bit array.
	 *
	 * @param i position of bit to set, relative to the entire
	 * bloom filter (not relative to the window!)
	 * @param val true or false
	 */
	void set(size_t i, bool val)
	{
		if (i >= m_startBitPos && i <= m_endBitPos)
			BloomFilter::set(i - m_startBitPos, val);
	}

	/** Add the object to this set. */
	void insert(const Bloom::key_type& key)
	{
		insert(Bloom::hash(key) % m_fullBloomSize);
	}

	/** Operator for reading a bloom filter from a stream. */
	friend std::istream& operator>>(std::istream& in, BloomFilterWindow& o)
	{
		o.read(in, Bloom::LOAD_OVERWRITE);
		return in;
	}

	/** Operator for writing the bloom filter to a stream. */
	friend std::ostream& operator<<(std::ostream& out, const BloomFilterWindow& o)
	{
		o.write(out);
		return out;
	}

	/** Read a bloom filter window from a stream. */
	void read(std::istream& in,
			Bloom::LoadType loadType = Bloom::LOAD_OVERWRITE,
			unsigned shrinkFactor = 1)
	{
		Bloom::FileHeader header = Bloom::readHeader(in);

		m_fullBloomSize = header.fullBloomSize;
		m_startBitPos = header.startBitPos;
		m_endBitPos = header.endBitPos;

		// alter the dimensions that we pass into Bloom::readData
		// so that we load the data into a bit array that is
		// exactly the size of the window (not the full size of the
		// containing bloom filter)

		header.fullBloomSize = header.endBitPos - header.startBitPos + 1;
		header.startBitPos = 0;
		header.endBitPos = header.fullBloomSize - 1;

		Bloom::readData(m_array, header, in, loadType, shrinkFactor);
	}

	/** Write a bloom filter window to a stream. */
	void write(std::ostream& out) const
	{
		Bloom::write(m_array, m_fullBloomSize, m_startBitPos,
			m_endBitPos, out);
	}

private:

	size_t m_fullBloomSize;
	size_t m_startBitPos, m_endBitPos;
};

#endif
