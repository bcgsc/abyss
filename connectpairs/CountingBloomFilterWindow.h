#ifndef COUNTINGBLOOMFILTERWINDOW_H
#define COUNTINGBLOOMFILTERWINDOW_H 1

#include "CountingBloomFilter.h"
#include <vector>

class CountingBloomFilterWindow : public CountingBloomFilter
{
  public:

	/** Constructor.
	 *
	 * @param bits size in bits of the containing counting bloom filter
	 * @param startBitPos index of first bit in the window
	 * @param endBitPos index of last bit in the window
	 */
	CountingBloomFilterWindow(size_t bits, size_t startBitPos, size_t endBitPos)
	{ 
		for (unsigned i = 0; i < MAX_COUNT; i++)
			m_data.push_back(new BloomFilterWindow(bits, startBitPos, endBitPos));
	}

};

#endif
