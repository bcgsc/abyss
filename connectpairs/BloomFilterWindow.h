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
	 * @param bits size in bits of the containing bloom filter
	 * @param startBitPos index of first bit in the window
	 * @param endBitPos index of last bit in the window
	 */
	BloomFilterWindow(size_t bits, size_t startBitPos, size_t endBitPos) :
		m_bits(bits),
		m_startBitPos(startBitPos),
		m_endBitPos(endBitPos)//,
	{
		assert(startBitPos <= endBitPos);
		m_array.resize(endBitPos - startBitPos + 1, 0);
	}

	/** Return whether the specified bit is set. */
	virtual bool operator[](size_t i) const
	{
		assert(i >= m_startBitPos && i <= m_endBitPos);
		return m_array[i - m_startBitPos];
	}

	/** Return whether the object is present in this set. */
	virtual bool operator[](const key_type& key) const
	{
		size_t i = hash(key) % m_bits;
		if (i < m_startBitPos || i > m_endBitPos)
			return false;
		return m_array[i - m_startBitPos];
	}

	/** Add the object with the specified index to this set. */
	virtual void insert(size_t i)
	{
		if (i < m_startBitPos || i > m_endBitPos)
			return;
		m_array[i - m_startBitPos] = true;
	}

protected:

	virtual void writeBloomDimensions(std::ostream& out) const
	{
		out << m_bits
			<< '\t' << m_startBitPos
			<< '\t' << m_endBitPos
			<< '\n';
	}

private:

	size_t m_bits;
	size_t m_startBitPos, m_endBitPos;
};

#endif
