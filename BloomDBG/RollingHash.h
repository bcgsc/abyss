#ifndef ABYSS_ROLLING_HASH_H
#define ABYSS_ROLLING_HASH_H 1

#include "lib/bloomfilter-2dfba08d120d7659e8c75cf5c501b3b9040e98cb/rolling.h"
#include <string>
#include <vector>
#include <cassert>
#include <boost/dynamic_bitset.hpp>

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
	 * a single seed hash value.
	 */
	std::vector<size_t> multiHash(size_t seedHash) const
	{
		std::vector<size_t> hashes(m_numHashes);
		for (unsigned i = 0; i < m_numHashes; i++) {
			hashes.at(i) = rol(varSeed, i) ^ seedHash;
		}
		return hashes;
	}

	/**
	 * Mask "don't care" positions in the current k-mer by
	 * replacing them with 'X' characters.
	 */
	void maskKmer()
	{
		assert(m_kmerMask.length() == m_k);
		/*
		 * note: dynamic_bitset stores the bits in reverse order relative
		 * to the string representation (e.g. "1101" is stored as "1011")
		 */
		for(size_t i = 0; i < m_kmerMask.length(); ++i) {
			if (m_kmerMask.at(i) == '0')
				m_kmer.at(i) = 'X';
		}
	}

	/**
	 * Restore the current k-mer to its unmasked state.
	 */
	void unmaskKmer()
	{
		m_kmer = m_unmaskedKmer;
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
		: m_numHashes(numHashes), m_k(k), m_hash1(0), m_rcHash1(0)
	{
		initKmerMask(std::string(k, '1'));
	}

	/**
	 * Constructor. Construct RollingHash object when initial k-mer
	 * is unknown.
	 * @param numHashes number of pseudo-independent hash values to compute
	 * for each k-mer
	 * @param k k-mer length
	 */
	RollingHash(unsigned numHashes, unsigned k, const std::string& kmerMask)
		: m_numHashes(numHashes), m_k(k), m_hash1(0), m_rcHash1(0)
	{
		initKmerMask(kmerMask);
	}

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
		initKmerMask(std::string(k, '1'));

		/* init rolling hash state */
		reset(kmer);
	}

	/**
	 * Constructor. Construct RollingHash object while specifying
	 * initial k-mer to be hashed.  This version of the constructor
	 * takes an additional bitmask argument that indicates which
	 * base positions should be ignored when hashing each k-mer.
	 * @param kmer initial k-mer for initializing hash value(s)
	 * @param numHashes number of pseudo-independent hash values to compute
	 * for each k-mer
	 * @param k k-mer length
	 * @param kmerMask bitmask indicating which base positions should
	 * be ignored during the hash calculation of each k-mer
	 */
	RollingHash(const std::string& kmer, unsigned numHashes, unsigned k,
		const std::string& kmerMask)
		: m_numHashes(numHashes), m_k(k), m_kmer(kmer), m_unmaskedKmer(kmer)
	{
		initKmerMask(kmerMask);

		/* init rolling hash state */
		reset(kmer);
	}

	/**
	 * Initialize k-mer mask (spaced seed).
	 */
	void initKmerMask(const std::string& kmerMask)
	{
		assert(kmerMask.length() == m_k);
		m_kmerMask = kmerMask;

		/* use the more efficient unmasked operations if k-mer mask is all "1"s */
		m_useKmerMask = false;
		for (size_t i = 0; i < m_kmerMask.length(); ++i) {
			if (m_kmerMask.at(i) == '0') {
				m_useKmerMask = true;
				break;
			}
		}
	}

	/**
	 * Initialize hash state from sequence.
	 * @param kmer k-mer used to initialize hash state
	 */
	void reset(const std::string& kmer)
	{
		if (m_useKmerMask)
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
		m_kmer = kmer;
		m_unmaskedKmer = kmer;

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
		m_hash1 = getFhval(m_kmer.c_str(), m_k);

		/* compute first hash function for reverse complement
		 * of k-mer */
		m_rcHash1 = getRhval(m_kmer.c_str(), m_k);

		/* compute hash values */
		m_hashes = multiHash(canonicalHash(m_hash1, m_rcHash1));

		/* restore k-mer to unmasked state */
		unmaskKmer();
	}

	/**
	 * Initialize hash values from sequence.
	 * @param kmer k-mer used to initialize hash state
	 */
	void resetUnmasked(const std::string& kmer)
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
	 * Compute hash values for next k-mer to the right and
	 * update internal state.
	 * @param charOut leftmost base of current k-mer
	 * @param charIn rightmost base of next k-mer
	 * @return vector of hash values for next k-mer
	 */
	void rollRight(unsigned char charOut, unsigned char charIn)
	{
		if (m_useKmerMask)
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
		std::rotate(m_kmer.begin(), m_kmer.begin() + 1, m_kmer.end());
		m_kmer.at(m_k - 1) = charIn;
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
		if (m_useKmerMask)
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
		std::rotate(m_kmer.rbegin(), m_kmer.rbegin() + 1, m_kmer.rend());
		m_kmer.at(0) = charIn;
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
		m_hashes = multiHash(seed);
	}

	/** Get hash values for current k-mer */
	const std::vector<size_t>& getHash() const
	{
		assert(!m_hashes.empty());
		return m_hashes;
	}

	/** Equality operator */
	bool operator==(const RollingHash& o) const
	{
		return
			m_numHashes == o.m_numHashes &&
			m_k == o.m_k &&
			m_hashes == o.m_hashes &&
			m_hash1 == o.m_hash1 &&
			m_rcHash1 == o.m_rcHash1;
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
	/** true if the k-mer mask should be used during hashing operations */
	bool m_useKmerMask;
	/** k-mer mask */
	std::string m_kmerMask;
	/** current k-mer (used only when k-mer mask is in effect) */
	std::string m_kmer;
	/** unmasked version of current k-mer (used only when k-mer mask is in effect) */
	std::string m_unmaskedKmer;

};

#endif
