#ifndef ABYSS_ROLLING_HASH_H
#define ABYSS_ROLLING_HASH_H 1

#include "config.h"
#include "lib/bloomfilter/rolling.h"
#include "BloomDBG/MaskedKmer.h"
#include <string>
#include <vector>
#include <cassert>
#include <boost/dynamic_bitset.hpp>
#include <cstring>

inline uint64_t spacedSeedHash(uint64_t& fhVal, uint64_t& rhVal,
	const char* kmerSeq, const char* spacedSeed, unsigned k)
{
	fhVal = 0;
	rhVal = 0;
    for(unsigned i=0; i<k; i++) {
		if (spacedSeed[i] != '0') {
			fhVal ^= rol(seedTab[(unsigned char)kmerSeq[i]], k-1-i);
			rhVal ^= rol(seedTab[(unsigned char)kmerSeq[i]+cpOff], i);
		}
	}
	return (rhVal < fhVal) ? rhVal : fhVal;
}

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

public:

	/**
	 * Default constructor.
	 */
	RollingHash() : m_k(0), m_hash1(0), m_rcHash1(0) {}

	/**
	 * Constructor. Construct RollingHash object when initial k-mer
	 * is unknown.
	 * @param numHashes number of pseudo-independent hash values to compute
	 * for each k-mer
	 * @param k k-mer length
	 */
	RollingHash(unsigned k) : m_k(k), m_hash1(0), m_rcHash1(0) {}

	/**
	 * Constructor. Construct RollingHash object while specifying
	 * initial k-mer to be hashed.
	 * @param kmer initial k-mer for initializing hash value(s)
	 * @param numHashes number of pseudo-independent hash values to compute
	 * for each k-mer
	 * @param k k-mer length
	 */
	RollingHash(const std::string& kmer, unsigned k)
		: m_k(k), m_hash1(0), m_rcHash1(0)
	{
		/* init rolling hash state */
		reset(kmer);
	}

	/**
	 * Initialize hash state from sequence.
	 * @param kmer k-mer used to initialize hash state
	 */
	void reset(const std::string& kmer)
	{
		if (!MaskedKmer::mask().empty())
			resetMasked(kmer);
		else
			resetUnmasked(kmer);
	}

	/**
	 * Initialize hash values from sequence. When computing the hash
	 * value, mask out "don't care" positions as per the active
	 * k-mer mask.
	 * @param kmer k-mer used to initialize hash state
	 */
	void resetMasked(const std::string& kmer)
	{
		assert(kmer.length() == m_k);

		/* store copy of k-mer for future rolling/masking ops */
		std::copy(kmer.c_str(), kmer.c_str() + m_k, m_kmer);

		resetMasked();
	}

	/**
	 * Initialize hash values from current k-mer. When computing the hash
	 * value, mask out "don't care" positions as per the active
	 * k-mer mask.
	 */
	void resetMasked()
	{
		spacedSeedHash(m_hash1, m_rcHash1, m_kmer,
			MaskedKmer::mask().c_str(), m_k);
	}

	/**
	 * Initialize hash values from sequence.
	 * @param kmer k-mer used to initialize hash state
	 */
	void resetUnmasked(const std::string& kmer)
	{
		/* compute first hash function for k-mer */
		m_hash1 = getFhval(kmer.c_str(), m_k);

		/* compute first hash function for reverse complement
		 * of k-mer */
		m_rcHash1 = getRhval(kmer.c_str(), m_k);
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
		if (!MaskedKmer::mask().empty())
			rollRightMasked(charOut, charIn);
		else
			rollRightUnmasked(charOut, charIn);
	}

	/**
	 * Compute hash values for next k-mer to the right and
	 * update internal state.  When computing the new hash, mask
	 * out "don't care" positions according to the active
	 * k-mer mask.
	 * @param charOut leftmost base of current k-mer
	 * @param charIn rightmost base of next k-mer
	 * @return vector of hash values for next k-mer
	 */
	void rollRightMasked(unsigned char, unsigned char charIn)
	{
		assert(m_k >= 2);
		memmove(m_kmer, m_kmer + 1, m_k - 1);
		m_kmer[m_k - 1] = charIn;
		resetMasked();
	}

	/**
	 * Compute hash values for next k-mer to the right and
	 * update internal state.
	 * @param charOut leftmost base of current k-mer
	 * @param charIn rightmost base of next k-mer
	 * @return vector of hash values for next k-mer
	 */
	void rollRightUnmasked(unsigned char charOut, unsigned char charIn)
	{
		/* update first hash function */
		rollHashesRight(m_hash1, m_rcHash1, charOut, charIn, m_k);
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
		if (!MaskedKmer::mask().empty())
			rollLeftMasked(charIn, charOut);
		else
			rollLeftUnmasked(charIn, charOut);
	}

	/**
	 * Compute hash values for next k-mer to the left and
	 * update internal state.  When computing the new hash, mask
	 * out "don't care" positions according to the active
	 * k-mer mask.
	 * @param charOut leftmost base of current k-mer
	 * @param charIn rightmost base of next k-mer
	 * @return vector of hash values for next k-mer
	 */
	void rollLeftMasked(unsigned char charIn, unsigned char)
	{
		assert(m_k >= 2);
		memmove(m_kmer + 1, m_kmer, m_k - 1);
		m_kmer[0] = charIn;
		resetMasked();
	}

	/**
	 * Compute hash values for next k-mer to the left and
	 * update internal state.
	 * @param charOut leftmost base of current k-mer
	 * @param charIn rightmost base of next k-mer
	 * @return vector of hash values for next k-mer
	 */
	void rollLeftUnmasked(unsigned char charIn, unsigned char charOut)
	{
		/* update first hash function */
		rollHashesLeft(m_hash1, m_rcHash1, charIn, charOut, m_k);
	}

	/** Get hash value for current k-mer */
	size_t getHash() const
	{
		return canonicalHash(m_hash1, m_rcHash1);
	}

	/** Equality operator */
	bool operator==(const RollingHash& o) const
	{
		return
			m_k == o.m_k &&
			m_hash1 == o.m_hash1 &&
			m_rcHash1 == o.m_rcHash1;
	}

private:

	/** k-mer length */
	unsigned m_k;
	/** value of first hash function for current k-mer */
	size_t m_hash1;
	/** value of first hash function for current k-mer, after
	 * reverse-complementing */
	size_t m_rcHash1;
	/** current k-mer (used only when k-mer mask is in effect) */
	char m_kmer[MAX_KMER];
};

#endif
