/**
 * A Bloom filter
 * Copyright 2013 Shaun Jackman
 */
#ifndef BLOOMFILTER_H
#define BLOOMFILTER_H 1

#include "Common/HashFunction.h"
#include "Common/Kmer.h"
#include <algorithm>
#include <vector>

/** A Bloom filter. */
class BloomFilter {
  public:
	/** The key type. */
	typedef Kmer key_type;

	/** Constructor. */
	BloomFilter(size_t n) : m_array(n) { }

	/** Return the hash value of this object. */
	static size_t hash(const key_type& key)
	{
		return hashmem(&key, sizeof key);
	}

	/** Return the size of the bit array. */
	size_t size() const { return m_array.size(); }

	/** Return the population count, i.e. the number of set bits. */
	size_t popcount() const
	{
		return std::count(m_array.begin(), m_array.end(), true);
	}

	/** Return whether the specified bit is set. */
	bool operator[](size_t i) const
	{
		assert(i < m_array.size());
		return m_array[i];
	}

	/** Return whether the object is present in this set. */
	bool operator[](const key_type& key) const
	{
		return m_array[hash(key) % m_array.size()];
	}

	/** Add the object with the specified index to this set. */
	void insert(size_t index)
	{
		assert(index < m_array.size());
		m_array[index] = true;
	}

	/** Add the object to this set. */
	void insert(const key_type& key)
	{
		m_array[hash(key) % m_array.size()] = true;
	}

  private:
	std::vector<bool> m_array;
};

#endif
