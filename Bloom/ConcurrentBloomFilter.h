#ifndef CONCURRENTBLOOMFILTER_H
#define CONCURRENTBLOOMFILTER_H

#ifndef _OPENMP
# error ConcurrentBloomFilter class requires a compiler that supports OpenMP
#endif

#include "config.h"
#include <vector>
#include <omp.h>

/**
 * A wrapper class that makes a Bloom filter
 * thread-safe.
 */
template <class BloomFilterType>
class ConcurrentBloomFilter
{

public:

	/** Constructor */
	ConcurrentBloomFilter(BloomFilterType& bloom, size_t numLocks,
		size_t hashSeed=0) : m_bloom(bloom), m_locks(numLocks),
		m_hashSeed(hashSeed)
	{
		m_windowSize = bloom.size() / numLocks;
		// round down to the nearest byte boundary,
		// because bytes that span locks will
		// cause concurrency issues
		m_windowSize -= m_windowSize % 8;
		assert(numLocks < bloom.size());
		for (size_t i = 0; i < m_locks.size(); i++)
			omp_init_lock(&(m_locks.at(i)));
	}

	/** Destructor */
	~ConcurrentBloomFilter()
	{
		for (size_t i = 0; i < m_locks.size(); i++)
			omp_destroy_lock(&(m_locks.at(i)));
	}

	/** Return whether the specified bit is set. */
	bool operator[](size_t i) const
	{
		assert(i < m_bloom.size());
		bool bit;
		getLock(i);
			bit = m_bloom[i];
		releaseLock(i);
		return bit;
	}

	/** Return whether the object is present in this set. */
	bool operator[](const Bloom::key_type& key) const
	{
		return *this[Bloom::hash(key, m_hashSeed) % m_bloom.size()];
	}

	/** Add the object with the specified index to this set. */
	void insert(size_t index)
	{
		assert(index < m_bloom.size());
		getLock(index);
			m_bloom.insert(index);
		releaseLock(index);
	}

	/** Add the object to this set. */
	void insert(const Bloom::key_type& key)
	{
		insert(Bloom::hash(key, m_hashSeed) % m_bloom.size());
	}

private:

	void getLock(size_t bitIndex)
	{
		assert(bitIndex < m_bloom.size());
		size_t lockIndex = std::min(bitIndex / m_windowSize, m_locks.size() - 1);
		omp_set_lock(&(m_locks.at(lockIndex)));
	}

	void releaseLock(size_t bitIndex)
	{
		assert(bitIndex < m_bloom.size());
		size_t lockIndex = std::min(bitIndex / m_windowSize, m_locks.size() - 1);
		omp_unset_lock(&(m_locks.at(lockIndex)));
	}

	BloomFilterType& m_bloom;
	std::vector<omp_lock_t> m_locks;
	size_t m_hashSeed;
	size_t m_windowSize;
};

#endif
