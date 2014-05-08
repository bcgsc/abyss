/**
 * A counting Bloom filter
 * Copyright 2014 bcgsc
 */
#ifndef COUNTINGBLOOMFILTER_H
#define COUNTINGBLOOMFILTER_H 1

#include "BloomFilter.h"
#include <vector>

/** A counting Bloom filter. */
template<class T>
class CountingBloomFilter : public BloomFilterBase
{
  public:

	/** The maximum count of an element in this multiset. */
	static const unsigned MAX_COUNT = 2;

	/** Constructor */
	CountingBloomFilter() {}

	/** Constructor */
	CountingBloomFilter(size_t n)
	{
		for (unsigned i = 0; i < MAX_COUNT; i++)
			m_data.push_back(new BloomFilter(n));
	}

	/** Destructor */
	virtual ~CountingBloomFilter()
	{
		typedef std::vector<BloomFilter*>::iterator Iterator;
		for (Iterator i = m_data.begin(); i != m_data.end(); i++) {
			assert(*i != NULL);
			delete *i;
		}
	}

	/** Return the size of the bit array. */
	size_t size() const
	{
		assert(m_data.back() != NULL);
		return m_data.back()->size();
	}

	/** Return the number of elements with count >= MAX_COUNT. */
	size_t popcount() const
	{
		assert(m_data.back() != NULL);
		return m_data.back()->popcount();
	}

	/** Return the estimated false positive rate */
	double FPR() const
	{
		return (double)popcount() / size();
	}

	/** Return whether the element with this index has count >=
	 * MAX_COUNT.
	 */
	bool operator[](size_t i) const
	{
		assert(m_data.back() != NULL);
		return (*m_data.back())[i];
	}

	/** Return whether this element has count >= MAX_COUNT. */
	bool operator[](const key_type& key) const
	{
		assert(m_data.back() != NULL);
		return (*m_data.back())[hash(key) % m_data.back()->size()];
	}

	/** Add the object with the specified index to this multiset. */
	virtual void insert(size_t index)
	{
		for (unsigned i = 0; i < MAX_COUNT; ++i) {
			assert(m_data.at(i) != NULL);
			if (!(*m_data[i])[index]) {
				m_data[i]->insert(index);
				break;
			}
		}
	}

	/** Add the object to this counting multiset. */
	virtual void insert(const key_type& key)
	{
		assert(m_data.back() != NULL);
		insert(hash(key) % m_data.back()->size());
	}

	/** Get the Bloom filter for a given level */
	BloomFilter& getBloomFilter(unsigned level)
	{
		assert(m_data.at(level) != NULL);
		return *m_data.at(level);
	}

	virtual void write(std::ostream& out) const
	{
		assert(m_data.back() != NULL);
		out << *m_data.back();
	}

  protected:
	std::vector<BloomFilter*> m_data;

};

#endif
