#ifndef ABYSS_ROLLING_HASH_H
#define ABYSS_ROLLING_HASH_H 1

#include "config.h"
#include "lib/rolling-hash/rolling.h"
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
		if (!MaskedKmer::mask().empty())
			resetMasked(kmer.c_str());
		else
			resetUnmasked(kmer);
	}

	/**
	 * Initialize hash values from current k-mer. When computing the hash
	 * value, mask out "don't care" positions as per the active
	 * k-mer mask.
	 */
	void resetMasked(const char* kmer)
	{
		const std::string& spacedSeed = MaskedKmer::mask();
		assert(spacedSeed.length() == m_k);

		/* compute first hash function for k-mer */
		uint64_t hash1 = getFhval(m_hash1, spacedSeed.c_str(), kmer, m_k);

		/* compute first hash function for reverse complement of k-mer */
		uint64_t rcHash1 = getRhval(m_rcHash1, spacedSeed.c_str(), kmer, m_k);

		m_hash = canonicalHash(hash1, rcHash1);
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

		m_hash = canonicalHash(m_hash1, m_rcHash1);
	}

	/**
	 * Compute hash values for next k-mer to the right and
	 * update internal state.
	 * @param kmer current k-mer
	 * @param nextKmer k-mer we are rolling into
	 */
	void rollRight(const char* kmer, char charIn)
	{
		if (!MaskedKmer::mask().empty())
			rollRightMasked(kmer, charIn);
		else
			rollRightUnmasked(kmer, charIn);
	}

	/**
	 * Compute hash values for next k-mer to the right and
	 * update internal state.  When computing the new hash, mask
	 * out "don't care" positions according to the active
	 * k-mer mask.
	 * @param kmer current k-mer
	 * @param nextKmer k-mer we are rolling into
	 */
	void rollRightMasked(const char* kmer, char charIn)
	{
		const std::string& spacedSeed = MaskedKmer::mask();
		m_hash = rollHashesRight(m_hash1, m_rcHash1, spacedSeed.c_str(),
			kmer, charIn, m_k);
	}

	/**
	 * Compute hash values for next k-mer to the right and
	 * update internal state.
	 * @param kmer current k-mer
	 * @param nextKmer k-mer we are rolling into
	 */
	void rollRightUnmasked(const char* kmer, char charIn)
	{
		/* update first hash function */
		rollHashesRight(m_hash1, m_rcHash1, kmer[0], charIn, m_k);
		m_hash = canonicalHash(m_hash1, m_rcHash1);
	}

	/**
	 * Compute hash values for next k-mer to the left and
	 * update internal state.
	 * @param prevKmer k-mer we are rolling into
	 * @param kmer current k-mer
	 */
	void rollLeft(char charIn, const char* kmer)
	{
		if (!MaskedKmer::mask().empty())
			rollLeftMasked(charIn, kmer);
		else
			rollLeftUnmasked(charIn, kmer);
	}

	/**
	 * Compute hash values for next k-mer to the left and
	 * update internal state.  When computing the new hash, mask
	 * out "don't care" positions according to the active
	 * k-mer mask.
	 * @param prevKmer k-mer we are rolling into
	 * @param kmer current k-mer
	 */
	void rollLeftMasked(char charIn, const char* kmer)
	{
		const std::string& spacedSeed = MaskedKmer::mask();
		m_hash = rollHashesLeft(m_hash1, m_rcHash1, spacedSeed.c_str(),
			kmer, charIn, m_k);
	}

	/**
	 * Compute hash values for next k-mer to the left and
	 * update internal state.
	 * @param prevKmer k-mer we are rolling into
	 * @param kmer current k-mer
	 */
	void rollLeftUnmasked(char charIn, const char* kmer)
	{
		/* update first hash function */
		rollHashesLeft(m_hash1, m_rcHash1, charIn, kmer[m_k-1], m_k);
		m_hash = canonicalHash(m_hash1, m_rcHash1);
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
		uint64_t tmpHashes[MAX_HASHES];
		multiHash(tmpHashes, m_hash, m_numHashes, m_k);
		for (unsigned i = 0; i < m_numHashes; ++i) {
			hashes[i] = (size_t)tmpHashes[i];
		}
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
	 * Set the base at a given position in the k-mer and update the hash
	 * value accordingly.
	 * @param kmer point to the k-mer char array
	 * @param pos position of the base to be changed
	 * @param base new value for the base
	 */
	void setBase(char* kmer, unsigned pos, char base)
	{
		if (!MaskedKmer::mask().empty())
			setBaseMasked(kmer, pos, base);
		else
			setBaseUnmasked(kmer, pos, base);
	}

	/**
	 * Set the base at a given position in the k-mer and update the hash
	 * value accordingly.
	 * @param kmer point to the k-mer char array
	 * @param pos position of the base to be changed
	 * @param base new value for the base
	 */
	void setBaseMasked(char* kmer, unsigned pos, char base)
	{
		const std::string& spacedSeed = MaskedKmer::mask();
		assert(spacedSeed.length() == m_k);
		m_hash = ::setBase(m_hash1, m_rcHash1, spacedSeed.c_str(), kmer,
			pos, base, m_k);
	}

	/**
	 * Set the base at a given position in the k-mer and update the hash
	 * value accordingly.
	 * @param kmer point to the k-mer char array
	 * @param pos position of the base to be changed
	 * @param base new value for the base
	 */
	void setBaseUnmasked(char* kmer, unsigned pos, char base)
	{
		m_hash = ::setBase(m_hash1, m_rcHash1, kmer, pos, base, m_k);
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
