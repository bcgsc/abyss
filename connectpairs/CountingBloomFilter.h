/**
 * A counting Bloom filter
 * Copyright 2013 Shaun Jackman
 */
#ifndef COUNTINGBLOOMFILTER_H
#define COUNTINGBLOOMFILTER_H 1

#include "BloomFilter.h"
#include <vector>

/** A counting Bloom filter. */
class CountingBloomFilter
{
  public:
	/** The key type. */
	typedef BloomFilter::key_type key_type;

	/** The maximum count of an element in this multiset. */
	static const unsigned MAX_COUNT = 2;

	/** Constructor */
	CountingBloomFilter(size_t n)
		: m_data(MAX_COUNT, BloomFilter(n)) { }

	/** Return the hash value of this object. */
	static size_t hash(const key_type& key)
	{
		return BloomFilter::hash(key);
	}

	/** Return the size of the bit array. */
	size_t size() const { return m_data.back().size(); }

	/** Return the number of elements with count >= MAX_COUNT. */
	size_t popcount() const
	{
		return m_data.back().popcount();
	}

	/** Return whether the element with this index has count >=
	 * MAX_COUNT.
	 */
	bool operator[](size_t i) const
	{
		assert(i < m_data.back().size());
		return m_data.back()[i];
	}

	/** Return whether this element has count >= MAX_COUNT. */
	bool operator[](const key_type& key) const
	{
		return m_data.back()[hash(key) % m_data.back().size()];
	}

	/** Add the object with the specified index to this multiset. */
	void insert(size_t index)
	{
		for (unsigned i = 0; i < MAX_COUNT; ++i) {
			if (!m_data[i][index]) {
				m_data[i].insert(index);
				break;
			}
		}
	}

	/** Add the object to this counting multiset. */
	void insert(const key_type& key)
	{
		insert(hash(key) % m_data.back().size());
	}

  private:
	std::vector<BloomFilter> m_data;
};

#endif
