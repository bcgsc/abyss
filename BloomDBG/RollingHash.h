#ifndef ABYSS_ROLLING_HASH_H
#define ABYSS_ROLLING_HASH_H 1

#include "config.h"

#include "BloomDBG/LightweightKmer.h"
#include "BloomDBG/MaskedKmer.h"
#include "Common/Sense.h"
#include "lib/nthash/nthash.hpp"

#include <algorithm>
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
	uint64_t canonicalHash(uint64_t hash, uint64_t rcHash) const
	{
		return (rcHash < hash) ? rcHash : hash;
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
	RollingHash(unsigned numHashes, unsigned k) : m_numHashes(numHashes),
		m_k(k), m_hash1(0), m_rcHash1(0) {}

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
		/* compute initial hash values for forward and reverse-complement k-mer */
		NTC64(kmer.c_str(), m_k, m_hash1, m_rcHash1);

		/* get canonical hash value from forward/reverse hash values */
		m_hash = canonicalHash(m_hash1, m_rcHash1);

		if (!MaskedKmer::mask().empty())
			m_hash = maskHash(m_hash1, m_rcHash1, MaskedKmer::mask().c_str(),
				kmer.c_str(), m_k);
	}

	/**
	 * Compute hash values for next k-mer to the right and
	 * update internal state.
	 * @param kmer current k-mer
	 * @param nextKmer k-mer we are rolling into
	 */
	void rollRight(const char* kmer, char charIn)
	{
		NTC64(kmer[0], charIn, m_k, m_hash1, m_rcHash1);
		m_hash = canonicalHash(m_hash1, m_rcHash1);

		if (!MaskedKmer::mask().empty()) {
			// TODO: copying the k-mer and shifting is very inefficient;
			// we need a specialized nthash function that rolls and masks
			// simultaneously
			LightweightKmer next(kmer);
			next.shift(SENSE, charIn);
			m_hash = maskHash(m_hash1, m_rcHash1, MaskedKmer::mask().c_str(),
				next.c_str(), m_k);
		}
	}

	/**
	 * Compute hash values for next k-mer to the left and
	 * update internal state.
	 * @param prevKmer k-mer we are rolling into
	 * @param kmer current k-mer
	 */
	void rollLeft(char charIn, const char* kmer)
	{
		NTC64L(kmer[m_k-1], charIn, m_k, m_hash1, m_rcHash1);
		m_hash = canonicalHash(m_hash1, m_rcHash1);

		if (!MaskedKmer::mask().empty()) {
			// TODO: copying the k-mer and shifting is very inefficient;
			// we need a specialized nthash function that rolls and masks
			// simultaneously
			LightweightKmer next(kmer);
			next.shift(ANTISENSE, charIn);
			m_hash = maskHash(m_hash1, m_rcHash1, MaskedKmer::mask().c_str(),
				next.c_str(), m_k);
		}
	}

	/**
	 * Get the seed hash value for the current k-mer. The seed hash
	 * value is used to calculate multiple pseudo-independant
	 * hash functions.
	 */
	size_t getHashSeed() const
	{
		return (size_t)m_hash;
	}

	/**
	 * Get hash values for current k-mer.
	 *
	 * @param hashes array for returned hash values
	 */
	void getHashes(size_t hashes[]) const
	{
		for (unsigned i = 0; i < m_numHashes; ++i)
			hashes[i] = NTE64(m_hash, m_k, i);
	}

	/** Equality operator */
	bool operator==(const RollingHash& o) const
	{
		/**
		 * Note: If hash seeds are equal, then the values
		 * for all hash functions will also be equal, since
		 * the hash values are calculated from the
		 * seed in a deterministic manner. In practice seed
		 * collision is very unlikely, though!
		 */
		return m_k == o.m_k && getHashSeed() == o.getHashSeed();
	}

	/** Inequality operator */
	bool operator!=(const RollingHash& o) const
	{
		return !(*this == o);
	}

	/**
	 * Change the hash value to reflect a change in the first/last base of
	 * the k-mer.
	 * @param kmer point to the k-mer char array
	 * @param dir if SENSE, change last base; if ANTISENSE,
	 * change first base
	 * @param base new value for the base
	 */
	void setLastBase(char* kmer, extDirection dir, char base)
	{
		if (dir == SENSE) {
			/* roll left to remove old last char */
			NTC64L(kmer[m_k-1], 'A', m_k, m_hash1, m_rcHash1);
			/* roll right to add new last char */
			NTC64('A', base, m_k, m_hash1, m_rcHash1);
		} else {
			/* roll right to remove old first char */
			NTC64(kmer[0], 'A', m_k, m_hash1, m_rcHash1);
			/* roll left to add new first char */
			NTC64L('A', base, m_k, m_hash1, m_rcHash1);
		}
		m_hash = canonicalHash(m_hash1, m_rcHash1);

		if (!MaskedKmer::mask().empty())
			m_hash = maskHash(m_hash1, m_rcHash1, MaskedKmer::mask().c_str(),
				kmer, m_k);
	}

private:

	/** number of hash functions */
	unsigned m_numHashes;
	/** k-mer length */
	unsigned m_k;
	/** value of first hash function for current k-mer */
	uint64_t m_hash1;
	/** value of first hash function for current k-mer, after
	 * reverse-complementing */
	uint64_t m_rcHash1;
	/** current canonical hash value */
	uint64_t m_hash;
};

#endif
