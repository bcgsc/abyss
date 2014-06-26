/**
 * A Bloom filter
 * Copyright 2013 Shaun Jackman
 */
#ifndef BLOOMFILTER_H
#define BLOOMFILTER_H 1

#include "Bloom/Bloom.h"
#include "Common/Kmer.h"
#include "Common/IOUtil.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <boost/dynamic_bitset.hpp>

/** A Bloom filter. */
class BloomFilter
{
  public:

	/** Constructor. */
	BloomFilter() { }

	/** Constructor. */
	BloomFilter(size_t n) : m_array(n) { }

	/** Return the size of the bit array. */
	size_t size() const { return m_array.size(); }

	/** Return the population count, i.e. the number of set bits. */
	size_t popcount() const
	{
		return m_array.count();
	}

	/** Return the estimated false positive rate */
	double FPR() const
	{
		return (double)popcount() / size();
	}

	/** Return whether the specified bit is set. */
	bool operator[](size_t i) const
	{
		assert(i < m_array.size());
		return m_array[i];
	}

	/** Return whether the object is present in this set. */
	bool operator[](const Bloom::key_type& key) const
	{
		return m_array[Bloom::hash(key) % m_array.size()];
	}

	/** Add the object with the specified index to this set. */
	void insert(size_t index)
	{
		set(index, true);
	}

	/** Set a bit to a specified value. */
	void set(size_t index, bool val)
	{
		assert(index < m_array.size());
		m_array[index] = val;
	}

	/** Add the object to this set. */
	void insert(const Bloom::key_type& key)
	{
		m_array[Bloom::hash(key) % m_array.size()] = true;
	}

	/** Operator for reading a bloom filter from a stream. */
	friend std::istream& operator>>(std::istream& in, BloomFilter& o)
	{
		o.read(in, Bloom::LOAD_OVERWRITE);
		return in;
	}

	/** Operator for writing the bloom filter to a stream. */
	friend std::ostream& operator<<(std::ostream& out, const BloomFilter& o)
	{
		o.write(out);
		return out;
	}

	/** Read a bloom filter from a stream. */
	void read(std::istream& in, Bloom::LoadType loadType = Bloom::LOAD_OVERWRITE,
		unsigned shrinkFactor = 1)
	{
		Bloom::read(m_array, in, loadType, shrinkFactor);
	}

	/** Write a bloom filter to a stream. */
	void write(std::ostream& out) const
	{
		Bloom::write(m_array, out);
	}

  protected:

	boost::dynamic_bitset<> m_array;
};

#endif
