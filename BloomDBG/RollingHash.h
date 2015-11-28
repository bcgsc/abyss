#ifndef ABYSS_ROLLING_HASH_H
#define ABYSS_ROLLING_HASH_H 1

#include "lib/bloomfilter-9061f087d8714109b865415f2ac05e4796d0cd74/rolling.h"
#include <string>
#include <vector>
#include <cassert>

class RollingHash
{
private:

	/**
	 * Determine the canonical hash value, given hash values for
	 * forward and reverse-complement of the same k-mer.
	 */
	size_t canonicalHash(size_t hash, size_t rcHash) const
	{
		return (rcHash < hash) ? rcHash : hash;
	}

	/** compute multiple pseudo-independent hash values using
	 * a single seed hash value */
	std::vector<size_t> multiHash(size_t seedHash) const
	{
		std::vector<size_t> hashes(m_numHashes);
		for (unsigned i = 0; i < m_numHashes; i++) {
			hashes.at(i) = rol(varSeed, i) ^ seedHash;
		}
		return hashes;
	}

public:

	/**
	 * Default constructor.
	 */
	RollingHash() : m_numHashes(0), m_k(0), m_hash1(0), m_rcHash1(0) {}

	/**
	 * Constructor. Construct RollingHash object when initial k-mer
	 * is unknown.
	 * @param numHashes number of pseudo-independent hash values to compute
	 * for each k-mer
	 * @param k k-mer length
	 */
	RollingHash(unsigned numHashes, unsigned k)
	: m_numHashes(numHashes), m_k(k), m_hash1(0), m_rcHash1(0) {}

	/**
	 * Constructor. Construct RollingHash object while specifying
	 * initial k-mer to be hashed.
	 * @param kmer initial k-mer for initializing hash value(s)
	 * @param numHashes number of pseudo-independent hash values to compute
	 * for each k-mer
	 * @param k k-mer length
	 */
	RollingHash(const std::string& kmer, unsigned numHashes, unsigned k)
		: m_numHashes(numHashes), m_k(k), m_hash1(0), m_rcHash1(0)
	{
		reset(kmer);
	}

	/** initialize hash values from seq */
	void reset(const std::string& kmer)
	{
		assert(kmer.length() == m_k);

		/* compute first hash function for k-mer */
		m_hash1 = getFhval(kmer.c_str(), m_k);

		/* compute first hash function for reverse complement
		 * of k-mer */
		m_rcHash1 = getRhval(kmer.c_str(), m_k);

		/* compute hash values */
		m_hashes = multiHash(canonicalHash(m_hash1, m_rcHash1));
	}

	/**
	 * Compute hash values for a neighbour k-mer on the right,
	 * without updating internal state.
	 * @param charOut leftmost base of current k-mer
	 * @param charIn rightmost base of k-mer to the right
	 * @return vector of hash values for next k-mer
	 */
	std::vector<size_t> peekRight(unsigned char charOut, unsigned char charIn) const
	{
		size_t hash1 = m_hash1;
		size_t rcHash1 = m_rcHash1;

		/* update first hash function */
		rollHashesRight(hash1, rcHash1, charOut, charIn, m_k);

		/* get seed value for computing rest of the hash functions */
		size_t seed = canonicalHash(hash1, rcHash1);

		/* compute hash values */
		return multiHash(seed);
	}

	/**
	 * Compute hash values for a neighbour k-mer on the left,
	 * without updating internal state.
	 * @param charIn leftmost base of k-mer to the left
	 * @param charOut rightmost base of current k-mer
	 * @return vector of hash values for next k-mer
	 */
	std::vector<size_t> peekLeft(unsigned char charIn, unsigned char charOut) const
	{
		size_t hash1 = m_hash1;
		size_t rcHash1 = m_rcHash1;

		/* update first hash function */
		rollHashesLeft(hash1, rcHash1, charIn, charOut, m_k);

		/* get seed value for computing rest of the hash functions */
		size_t seed = canonicalHash(hash1, rcHash1);

		/* compute hash values */
		return multiHash(seed);
	}

	/**
	 * Compute hash values for next k-mer to the right and
	 * update internal state.
	 * @param charOut leftmost base of current k-mer
	 * @param charIn rightmost base of next k-mer
	 * @return vector of hash values for next k-mer
	 */
	void rollRight(unsigned char charOut, unsigned char charIn)
	{
		/* update first hash function */
		rollHashesRight(m_hash1, m_rcHash1, charOut, charIn, m_k);

		/* get seed value for computing rest of the hash functions */
		size_t seed = canonicalHash(m_hash1, m_rcHash1);

		/* compute hash values */
		m_hashes = multiHash(seed);
	}

	/**
	 * Compute hash values for next k-mer to the left and
	 * update internal state.
	 * @param charOut leftmost base of current k-mer
	 * @param charIn rightmost base of next k-mer
	 * @return vector of hash values for next k-mer
	 */
	void rollLeft(unsigned char charIn, unsigned char charOut)
	{
		/* update first hash function */
		rollHashesLeft(m_hash1, m_rcHash1, charIn, charOut, m_k);

		/* get seed value for computing rest of the hash functions */
		size_t seed = canonicalHash(m_hash1, m_rcHash1);

		/* compute hash values */
		m_hashes = multiHash(seed);
	}

	/** Get hash values for current k-mer */
	const std::vector<size_t>& getHash() const
	{
		assert(!m_hashes.empty());
		return m_hashes;
	}

private:

	/** number of hash functions to compute at each position */
	unsigned m_numHashes;
	/** k-mer length */
	unsigned m_k;
	/** current values for each hash function */
	std::vector<size_t> m_hashes;
	/** value of first hash function for current k-mer */
	size_t m_hash1;
	/** value of first hash function for current k-mer, after
	 * reverse-complementing */
	size_t m_rcHash1;
};

#endif
