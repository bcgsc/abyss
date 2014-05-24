/**
 * A counting Bloom filter
 * Copyright 2014 bcgsc
 */
#ifndef COUNTINGBLOOMFILTER_H
#define COUNTINGBLOOMFILTER_H 1

#include "Common/BloomFilterBase.h"
#include <vector>
#include <math.h>

/** A counting Bloom filter. */
template<typename NumericType>
class CountingBloomFilter : public BloomFilterBase
{
  public:

	/** Constructor */
	CountingBloomFilter() {}

	/** Constructor */
	CountingBloomFilter(size_t n)
	{
		m_data = new std::vector<NumericType>(n);
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

	/** Return the number of elements with count >= MAX_COUNT. */
	size_t popcount() const
	{
		return uniqueEntries;
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
	NumericType operator[](size_t i) const
	{
		return m_data[i];
	}

	/** Return the count of this element. */
	NumericType operator[](const key_type& key) const
	{
		NumericType currentMin = m_data[hash(key, 0) % m_data.size()];
		for (unsigned int i = 1; i < hashNum; ++i) {
			NumericType min = m_data[hash(key, i) % m_data.size()];
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
		NumericType minEle = this[key];

		//update only those elements
		NumericType currentMin = m_data[hash(key, 0) % m_data.size()];
		for (unsigned int i = 1; i < hashNum; ++i) {
			size_t hashVal = hash(key, i) % m_data.size();
			NumericType val = m_data[hashVal];
			if (minEle == val) {
				insert(hashVal);
			}
		}
		if(minEle == 0)
			++uniqueEntries;
		else
			++replicateEntries;
	}

	virtual void write(std::ostream& out) const
	{
		assert(m_data != NULL);
		out << *m_data;
	}

	//overloaded from BloomFilterBase
	//need to impliment tracking of directionality
//	void loadSeq(unsigned k, const std::string& seq)
//	{
//		if (seq.size() < k)
//			return;
//		for (size_t i = 0; i < seq.size() - k + 1; ++i) {
//			std::string kmer = seq.substr(i, k);
//			size_t pos = kmer.find_last_not_of("ACGTacgt");
//			if (pos == std::string::npos) {
//				insert(Kmer(kmer));
//			} else
//				i += pos;
//		}
//	}

  protected:
	std::vector<T> m_data;
	unsigned hashNum;
	size_t uniqueEntries;
	size_t replicateEntries;

};

#endif
