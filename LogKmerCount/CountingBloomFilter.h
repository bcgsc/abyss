/**
 * A counting Bloom filter
 * Copyright 2014 bcgsc
 */
#ifndef COUNTINGBLOOMFILTER_H
#define COUNTINGBLOOMFILTER_H 1

#include "BloomFilter.h"
#include <vector>
#include <math.h>

/** A counting Bloom filter. */
template<class T>
class CountingBloomFilter : public BloomFilterBase
{
  public:

//	/** The maximum count of an element in this multiset. */
//	static const unsigned MAX_COUNT = 2;

	/** Constructor */
	CountingBloomFilter() {}

	/** Constructor */
	CountingBloomFilter(size_t n)
	{
		m_data = new std::vector<T>(n);
	}

	/** Destructor */
	virtual ~CountingBloomFilter()
	{
		delete m_data;
	}

	/** Return the size (in discrete elements) of the bit array. */
	size_t size() const
	{
		return m_data.size();
	}

	/** Return the estimated false positive rate */
	double FPR() const
	{
		return pow(1.0 - pow(1.0 - 1.0 / double(m_data.size()),
				double(uniqueEntries) * hashNum),
				double(hashNum));
	}

	/** Return the count of the single element (debugging purposes)
	 */
	bool operator[](size_t i) const
	{
		return m_data[i];
	}

	/** Return the count of this element. */
	T operator[](const key_type& key) const
	{
		T currentMin = m_data[hash(key, 0) % m_data.size()];
		for (unsigned int i = 1; i < hashNum; ++i) {
			T min = m_data[hash(key, i) % m_data.size()];
			if (min < currentMin) {
				currentMin = min;
			}
			if (0 == currentMin) {
				break;
			}
		}
		return currentMin;
	}

	/** Add the object with the specified index (debugging purposes). */
	virtual void insert(size_t index)
	{
		++m_data[index];
	}

	/** Add the object to this counting multiset.
	 *  If all values are the same update all
	 *  If some values are larger only update smallest counts*/
	virtual void insert(const key_type& key)
	{

		//check for which elements to update
		for (unsigned int i = 1; i < hashNum; ++i)
		{
			if()
			insert(hash(key, i) % m_data.size());
		}

		//update only those elements

		T currentMin = m_data[hash(key, 0) % m_data.size()];
		for (unsigned int i = 1; i < hashNum; ++i)
		{
			if()
			insert(hash(key, i) % m_data.size());
		}
	}

	virtual void write(std::ostream& out) const
	{
		assert(m_data != NULL);
		out << *m_data;
	}

  protected:
	std::vector<T> m_data;
	unsigned hashNum;
	unsigned uniqueEntries;

};

#endif
