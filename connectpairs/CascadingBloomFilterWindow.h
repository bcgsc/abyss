#ifndef COUNTINGBLOOMFILTERWINDOW_H
#define COUNTINGBLOOMFILTERWINDOW_H 1

#include "connectpairs/CascadingBloomFilter.h"
#include <vector>

class CascadingBloomFilterWindow : public CascadingBloomFilter
{
  public:

	/** Constructor.
	 *
	 * @param fullBloomSize size in bits of the containing counting bloom filter
	 * @param startBitPos index of first bit in the window
	 * @param endBitPos index of last bit in the window
	 */
	CascadingBloomFilterWindow(size_t fullBloomSize, size_t startBitPos, size_t endBitPos)
		: m_fullBloomSize(fullBloomSize)
	{
		for (unsigned i = 0; i < MAX_COUNT; i++)
			m_data.push_back(new BloomFilterWindow(fullBloomSize, startBitPos, endBitPos));
	}

	/** Add the object with the specified index to this multiset. */
	virtual void insert(size_t i)
	{
		CountingBloomFilter::insert(i);
	}

	/** Add the object to this counting multiset. */
	void insert(const key_type& key)
	{
		insert(hash(key) % m_fullBloomSize);
	}

  private:

	size_t m_fullBloomSize;
};

#endif
