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

	/**
	 * Compute multiple pseudo-independent hash values using
	 * a single seed hash value
	 */
	void multiHash(size_t seedHash)
	{
		for (unsigned i = 0; i < m_numHashes; i++) {
#if 1
			m_hashes[i] = rol(varSeed, i) ^ seedHash;
#else
			/* Hamid's version (has compile errors) */
			m_hashes[i] = seedHash * (i ^ m_kmer * varSeed);
			m_hashes[i] ^= m_hashes[i] >> varShift;
#endif
		}
	}

	/**
	 * Mask "don't care" positions in the current k-mer by
	 * replacing them with 'X' characters.
	 */
	void maskKmer()
	{
		assert(MaskedKmer::mask().length() == m_k);
		for(size_t i = 0; i < m_k; ++i) {
			if (MaskedKmer::mask().at(i) == '0')
				m_kmer[i] = 'X';
		}
	}

	/**
	 * Restore the current k-mer to its unmasked state.
	 */
	void unmaskKmer()
	{
		std::copy(m_unmaskedKmer, m_unmaskedKmer + m_k, m_kmer);
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
		std::copy(kmer.c_str(), kmer.c_str() + m_k, m_unmaskedKmer);

		resetMasked();
	}

	/**
	 * Initialize hash values from current k-mer. When computing the hash
	 * value, mask out "don't care" positions as per the active
	 * k-mer mask.
	 */
	void resetMasked()
	{
		/* replace "don't care" positions with 'X' */
		maskKmer();

		/* compute first hash function for k-mer */
		m_hash1 = getFhval(m_kmer, m_k);

		/* compute first hash function for reverse complement
		 * of k-mer */
		m_rcHash1 = getRhval(m_kmer, m_k);

		/* compute hash values */
		multiHash(canonicalHash(m_hash1, m_rcHash1));

		/* restore k-mer to unmasked state */
		unmaskKmer();
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

		/* compute hash values */
		multiHash(canonicalHash(m_hash1, m_rcHash1));
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
		memmove(m_unmaskedKmer, m_unmaskedKmer + 1, m_k - 1);
		m_unmaskedKmer[m_k - 1] = charIn;
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

		/* get seed value for computing rest of the hash functions */
		size_t seed = canonicalHash(m_hash1, m_rcHash1);

		/* compute hash values */
		multiHash(seed);
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
		memmove(m_unmaskedKmer + 1, m_unmaskedKmer, m_k - 1);
		m_unmaskedKmer[0] = charIn;
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

		/* get seed value for computing rest of the hash functions */
		size_t seed = canonicalHash(m_hash1, m_rcHash1);

		/* compute hash values */
		multiHash(seed);
	}

	/** Get hash values for current k-mer */
	const size_t* getHash() const
	{
		return m_hashes;
	}

	/** Equality operator */
	bool operator==(const RollingHash& o) const
	{
		return
			m_numHashes == o.m_numHashes &&
			m_k == o.m_k &&
			std::equal(m_hashes, m_hashes + m_numHashes, o.m_hashes) &&
			m_hash1 == o.m_hash1 &&
			m_rcHash1 == o.m_rcHash1;
	}

private:

	/** number of hash functions to compute at each position */
	unsigned m_numHashes;
	/** current values for each hash function */
	size_t m_hashes[MAX_HASHES];
	/** k-mer length */
	unsigned m_k;
	/** value of first hash function for current k-mer */
	size_t m_hash1;
	/** value of first hash function for current k-mer, after
	 * reverse-complementing */
	size_t m_rcHash1;
	/** current k-mer (used only when k-mer mask is in effect) */
	char m_kmer[MAX_KMER];
	/** unmasked version of current k-mer (used only when k-mer mask is in effect) */
	char m_unmaskedKmer[MAX_KMER];
};

#endif
