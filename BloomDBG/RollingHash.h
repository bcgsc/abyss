#ifndef ABYSS_ROLLING_HASH_H
#define ABYSS_ROLLING_HASH_H 1

#include "lib/bloomfilter-521e80c5c619a9a8e3d6389dc3b597a75bdf2aaa/rolling.h"
#include <string>
#include <vector>
#include <cassert>

class RollingHash
{
private:

	/**
	 * Initialize forward/reverse hash values for hash function 1.
	 * @param kmer sequence to be hashed
	 */
	void resetHash1(const std::string& kmer)
	{
		/* compute first hash function for k-mer */
		m_hash1 = getFhval(kmer.c_str(), m_k);

		/* compute first hash function for reverse complement
		 * of k-mer */
		m_rcHash1 = getRhval(kmer.c_str(), m_k);
	}

	/**
	 * Compute first hash value for next k-mer, for both
	 * possible orientations.
	 * @param charOut leftmost base of current k-mer
	 * @param charIn rightmost base of next k-mer
	 * @return vector of hash values for next k-mer
	 */
	void rollHash1(unsigned char charOut, unsigned char charIn)
	{
		/* roll first hash value for forward k-mer */
		m_hash1 = rol(m_hash1, 1)
			^ rol(seedTab[charOut], m_k)
			^ seedTab[charIn];

		/* roll first hash value for reverse complement k-mer */
		m_rcHash1 = ror(m_rcHash1, 1)
			^ ror(seedTab[charOut + cpOff], 1)
			^ rol(seedTab[charIn + cpOff], m_k - 1);
	}

	/** determine canonical value for first hash function */
	size_t canonicalHash1()
	{
		return (m_hash1 < m_rcHash1) ? m_hash1 : m_rcHash1;
	}

	/** compute multiple pseudo-independent hash values using
	 * a single seed hash value */
	std::vector<size_t> multiHash(size_t seedHash)
	{
		std::vector<size_t> hashes(m_numHashes);
		for (unsigned i = 0; i < m_numHashes; i++) {
			hashes.at(i) = rol(varSeed, i) ^ seedHash;
		}
		return hashes;
	}

public:

	/**
	 * Constructor.
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

		/* compute forward/reverse values for first hash function */
		resetHash1(kmer);

		/* compute hash values */
		m_hashes = multiHash(canonicalHash1());
	}

	/**
	 * Compute hash values for next k-mer without updating internal state.
	 * @param charOut leftmost base of current k-mer
	 * @param charIn rightmost base of next k-mer
	 * @return vector of hash values for next k-mer
	 */
	std::vector<size_t> peekRight(unsigned char charOut, unsigned char charIn)
	{
		/* save state */
		size_t hash1Orig = m_hash1;
		size_t rcHash1Orig = m_rcHash1;

		/* update first hash function */
		rollHash1(charOut, charIn);

		/* get seed value for computing rest of the hash functions */
	    size_t seed = canonicalHash1();

		/* restore state */
		m_hash1 = hash1Orig;
		m_rcHash1 = rcHash1Orig;

		/* compute hash values */
		return multiHash(seed);
	}

	/**
	 * Compute hash values for next k-mer and update internal state.
	 * @param charOut leftmost base of current k-mer
	 * @param charIn rightmost base of next k-mer
	 * @return vector of hash values for next k-mer
	 */
	void rollRight(unsigned char charOut, unsigned char charIn)
	{
		/* update first hash function */
		rollHash1(charOut, charIn);

		/* compute hash values */
		m_hashes = multiHash(canonicalHash1());
	}

	/** Get hash values for current k-mer */
	std::vector<size_t> getHash()
	{
		return m_hashes;
	}

private:

	/** number of hash functions to compute at each position */
	const unsigned m_numHashes;
	/** k-mer length */
	const unsigned m_k;
	/** current values for each hash function */
	std::vector<size_t> m_hashes;
	/** value of first hash function for current k-mer */
	size_t m_hash1;
	/** value of first hash function for current k-mer, after
	 * reverse-complementing */
	size_t m_rcHash1;
};

#endif
