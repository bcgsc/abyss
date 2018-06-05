/**
 * A cascading Bloom filter
 * Copyright 2013 Shaun Jackman
 */
#ifndef CascadingBLOOMFILTER_H
#define CascadingBLOOMFILTER_H 1

#include "Bloom/Bloom.h"
#include "BloomFilter.h"
#include <vector>

/** A Cascading Bloom filter. */
class CascadingBloomFilter
{
  public:

	/** Constructor */
	CascadingBloomFilter() {}

	/** Constructor */
	CascadingBloomFilter(size_t n, size_t max_count, size_t hashSeed=0) : m_hashSeed(hashSeed)
	{
		m_data.reserve(max_count);
		for (unsigned i = 0; i < max_count; i++)
			m_data.push_back(new Konnector::BloomFilter(n, hashSeed));
	}

	/** Destructor */
	~CascadingBloomFilter()
	{
		typedef std::vector<Konnector::BloomFilter*>::iterator Iterator;
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

	/** Return the number of elements with count >= max_count. */
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
	 * max_count.
	 */
	bool operator[](size_t i) const
	{
		assert(m_data.back() != NULL);
		return (*m_data.back())[i];
	}

	/** Return whether this element has count >= max_count. */
	bool operator[](const Bloom::key_type& key) const
	{
		assert(m_data.back() != NULL);
		return (*m_data.back())[Bloom::hash(key, m_hashSeed) % m_data.back()->size()];
	}

	/** Add the object with the specified index to this multiset. */
	void insert(size_t index)
	{
		for (unsigned i = 0; i < m_data.size(); ++i) {
			assert(m_data.at(i) != NULL);
			if (!(*m_data[i])[index]) {
				m_data[i]->insert(index);
				break;
			}
		}
	}

	/** Add the object to this Cascading multiset. */
	void insert(const Bloom::key_type& key)
	{
		assert(m_data.back() != NULL);
		insert(Bloom::hash(key, m_hashSeed) % m_data.back()->size());
	}

	/** Get the Bloom filter for a given level */
	Konnector::BloomFilter& getBloomFilter(unsigned level)
	{
		assert(m_data.at(level) != NULL);
		return *m_data.at(level);
	}

	void write(std::ostream& out) const
	{
		assert(m_data.back() != NULL);
		out << *m_data.back();
	}

	/** Operator for writing the bloom filter to a stream */
	friend std::ostream& operator<<(std::ostream& out, const CascadingBloomFilter& o)
	{
		o.write(out);
		return out;
	}

  private:
	size_t m_hashSeed;
	std::vector<Konnector::BloomFilter*> m_data;

};

#endif
