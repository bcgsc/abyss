/**
 * A cascading Bloom filter
 * Copyright 2015 Shaun Jackman, Ben Vandervalk.
 */
#ifndef HASH_AGNOSTIC_CASCADING_BLOOM_H
#define HASH_AGNOSTIC_CASCADING_BLOOM_H 1

#include "vendor/btl_bloomfilter/BloomFilter.hpp"
#include <vector>

/**
 * An implementation of a Cascading Bloom filter.
 * A Cascading Bloom filter implements a crude
 * counting mechanism using an array of _l_ Bloom
 * filters; we say that such a Bloom filter has
 * l _levels_. Each time an element is inserted, we
 * check for its presence in each level, and then
 * insert the element into the first Bloom filter
 * where the element is not already present.
 *
 * We use the Cascading Bloom filter to filter
 * out error k-mers from the de Bruijn graph, since
 * these k-mers typically only occur once in the
 * the data.
 */
class HashAgnosticCascadingBloom
{
  public:
	typedef uint64_t hash_t;
	/** Default constructor */
	HashAgnosticCascadingBloom()
	  : m_k(0)
	  , m_hashes(0)
	{}

	/**
	 * Constructor.
	 * @param size size of the Bloom filters (in bits)
	 * @param hashes number of hash functions
	 * @param levels number of levels in Cascading Bloom filter
	 * @param k k-mer size
	 */
	HashAgnosticCascadingBloom(size_t size, unsigned hashes, size_t levels, unsigned k)
	  : m_k(k)
	  , m_hashes(hashes)
	{
		m_data.reserve(levels);
		for (unsigned i = 0; i < levels; i++)
			m_data.push_back(new BloomFilter(size, hashes, k));
	}

	/**
	 * Constructor to load a single-level BloomFilter from
	 * files.  This is used to make BloomFilter support the
	 * same interface as HashAgnosticCascadingBloom.
	 */
	HashAgnosticCascadingBloom(const string& bloomPath) { loadFilter(bloomPath); }

	/** Destructor */
	~HashAgnosticCascadingBloom() { clear(); }

	/** Return k-mer size used by Bloom filter. */
	unsigned getKmerSize() const { return m_k; }

	/** Return number of hash functions used by Bloom filter */
	unsigned getHashNum() const { return m_hashes; }

	/** Return the size of the bit array. */
	size_t size() const
	{
		assert(m_data.back() != NULL);
		return m_data.back()->getFilterSize();
	}

	/** Return the size of the bit array. */
	size_t getFilterSize() const { return size(); }

	/** Return the number of elements with count >= levels. */
	size_t popcount() const
	{
		assert(m_data.back() != NULL);
		return m_data.back()->getPop();
	}

	/** Return number of levels in cascading Bloom filter */
	unsigned levels() const { return m_data.size(); }

	/** Return the estimated false positive rate */
	double FPR() const { return pow((double)popcount() / size(), m_hashes); }

	/**
	 * Return true if the element with the given hash values
	 * has count >= levels.
	 */
	bool contains(const std::vector<hash_t>& hashes) const
	{
		assert(m_data.back() != NULL);
		return m_data.back()->contains(hashes);
	}

	/**
	 * Return true if the element with the given hash values
	 * has count >= levels.
	 */
	bool contains(const hash_t hashes[]) const
	{
		assert(m_data.back() != NULL);
		return m_data.back()->contains(hashes);
	}

	/** Add the object with the specified index to this multiset. */
	void insert(const std::vector<hash_t>& hashes)
	{
		for (unsigned i = 0; i < m_data.size(); ++i) {
			assert(m_data.at(i) != NULL);
			if (!(*m_data[i]).contains(hashes)) {
				m_data[i]->insert(hashes);
				break;
			}
		}
	}

	/** Add the object with the specified index to this multiset. */
	void insert(const hash_t hashes[])
	{
		for (unsigned i = 0; i < m_data.size(); ++i) {
			assert(m_data.at(i) != NULL);
			if (!(*m_data[i]).contains(hashes)) {
				m_data[i]->insert(hashes);
				break;
			}
		}
	}

	/** Get the Bloom filter for a given level */
	BloomFilter& getBloomFilter(unsigned level)
	{
		assert(m_data.at(level) != NULL);
		return *m_data.at(level);
	}

	/** Operator for writing the Bloom filter to a stream */
	friend std::ostream& operator<<(std::ostream& out, const HashAgnosticCascadingBloom& o)
	{
		assert(o.m_data.size() > 0);
		assert(o.m_data.back() != NULL);
		/* o.m_data.back()->storeFilter(out); */
		out << *o.m_data.back();
		return out;
	}

	/** Load a Bloom filter from a file */
	void loadFilter(const string& bloomPath)
	{
		clear();
		BloomFilter* bloom = new BloomFilter(bloomPath);
		m_k = bloom->getKmerSize();
		m_hashes = bloom->getHashNum();
		m_data.push_back(bloom);
	}

  private:
	/** Free all allocated memory and reset parameters to defaults */
	void clear()
	{
		m_k = 0;
		m_hashes = 0;
		typedef std::vector<BloomFilter*>::iterator Iterator;
		for (Iterator i = m_data.begin(); i != m_data.end(); i++) {
			assert(*i != NULL);
			delete *i;
		}
		m_data.clear();
	}

	/** k-mer length */
	unsigned m_k;
	/** number of hash functions */
	unsigned m_hashes;
	/** the array of Bloom filters */
	std::vector<BloomFilter*> m_data;
};

#endif
