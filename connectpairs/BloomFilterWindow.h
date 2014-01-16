#ifndef BLOOMFILTERWINDOW_H
#define BLOOMFILTERWINDOW_H 1

#include "BloomFilterBase.h"
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

	/** Return whether the specified bit is set. */
	virtual bool operator[](size_t i) const
	{
		if (i >= m_startBitPos && i <= m_endBitPos)
			return BloomFilter::operator[](i - m_startBitPos);
		return false;
	}

	/** Return whether the object is present in this set. */
	virtual bool operator[](const key_type& key) const
	{
		return (*this)[hash(key) % m_fullBloomSize];
	}

	/** Add the object with the specified index to this set. */
	virtual void insert(size_t i)
	{
		if (i >= m_startBitPos && i <= m_endBitPos)
			BloomFilter::insert(i - m_startBitPos);
	}

	/** Add the object to this set. */
	virtual void insert(const key_type& key)
	{
		insert(hash(key) % m_fullBloomSize);
	}

	/**
	 * Return the full size of the containing bloom
	 * filter (in bits).
	 */
	size_t getFullBloomSize() const
	{
		return m_fullBloomSize;
	}

protected:

	virtual void writeBloomDimensions(std::ostream& out) const
	{
		out << m_fullBloomSize
			<< '\t' << m_startBitPos
			<< '\t' << m_endBitPos
			<< '\n';
	}

private:

	size_t m_fullBloomSize;
	size_t m_startBitPos, m_endBitPos;
};

#endif
