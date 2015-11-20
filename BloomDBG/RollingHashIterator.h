#ifndef ROLLING_HASH_ITERATOR_H
#define ROLLING_HASH_ITERATOR_H 1

#include <cstring>
#include <vector>
#include <cassert>
#include <limits>
#include "lib/bloomfilter-521e80c5c619a9a8e3d6389dc3b597a75bdf2aaa/rolling.h"

/**
 * Permitted characters in k-mers. All k-mers containing
 * other characters will be skipped.
 */
#define ACGT_CHARS "acgtACGT"

/**
 * Iterate over hash values for k-mers in a
 * given DNA sequence.
 *
 * This implementation uses a rolling hash
 * function to efficiently calculate
 * hash values for successive k-mers.
 */
class RollingHashIterator
{
public:

	/** data type for hash values */
	typedef uint64_t hash_t;

private:

	/**
	 * Init hash values from a k-mer.
	 */
	void initHash(const char* kmer)
	{
		/* compute first hash value for k-mer */
		m_hash1 = getFhval(kmer, m_k);

		/* compute first hash value for reverse complement
		 * of k-mer */
		m_rcHash1 = getRhval(kmer, m_k);

		/* determine "canonical" k-mer orientation */
		hash_t hash = (m_hash1 < m_rcHash1) ? m_hash1 : m_rcHash1;

		/* compute remaining hash values for k-mer */
		m_hashes.clear();
		for (size_t i = 0; i < m_numHashes; ++i)
			m_hashes.push_back(rol(varSeed, i) ^ hash);
	}

	/**
	 * Compute hash values for current k-mer using the
	 * hash values of k-mer immediately to the left.
	 */
	void rollHash()
	{
		assert(m_hashes.size() == m_numHashes);
		assert(m_pos > 0);
		assert(m_pos + m_k <= m_seqLen);

		const unsigned char charOut = m_seq[m_pos - 1];
		const unsigned char charIn = m_seq[m_pos + m_k - 1];

		/* roll first hash value for forward k-mer */
		m_hash1 = rol(m_hash1, 1)
			^ rol(seedTab[charOut], m_k)
			^ seedTab[charIn];

		/* roll first hash value for reverse complement k-mer */
		m_rcHash1 = ror(m_rcHash1, 1)
			^ ror(seedTab[charOut + cpOff], 1)
			^ rol(seedTab[charIn + cpOff], m_k - 1);

		/* determine canonical hash value */
		uint64_t hash = (m_hash1 < m_rcHash1) ? m_hash1 : m_rcHash1;

		/* compute remaining hash values */
		for (unsigned i = 0; i < m_numHashes; i++) {
			m_hashes.at(i) = rol(varSeed, i) ^ hash;
		}
	}

	/**
	 * Advance iterator right to the next valid k-mer.
	 */
	void next()
	{
		if (m_seqLen < m_k) {
			m_pos = std::numeric_limits<std::size_t>::max();
			return;
		}

		while(m_pos < m_seqLen - m_k + 1) {
			/* skip over k-mers with non-ACGT chars */
			if (m_nextInvalidChar - m_pos < m_k) {
				m_pos = m_nextInvalidChar + 1;
				m_nextInvalidChar = m_pos + strspn(m_seq + m_pos, ACGT_CHARS);
				m_hashes.clear();
			} else {
				/* we are positioned at the next valid k-mer */
				if (m_hashes.empty()) {
					/* we don't have hash values for the
					 * preceding k-mer, so we must compute
					 * the hash values from scratch */
					initHash(m_seq + m_pos);
				} else {
					/* compute new hash values based on
					 * hash values of preceding k-mer */
					assert(m_pos > 0);
					rollHash();
				}
				return;
			}
		}
		/* there are no more valid k-mers */
		m_pos = std::numeric_limits<std::size_t>::max();
	}

public:

	/**
	 * Default constructor. Creates an iterator pointing to
	 * the end of the iterator range.
	 */
	RollingHashIterator() : m_seq(NULL), m_seqLen(0), m_k(0),
		m_numHashes(0), m_pos(std::numeric_limits<std::size_t>::max()),
		m_hash1(0), m_rcHash1(0) {}

	/**
	 * Constructor.
	 * @param seq DNA sequence to be hashed
	 * @param k k-mer size
	 * @param numHashes number of hash values to compute
	 * for each k-mer
	 */
	RollingHashIterator(const char* seq, unsigned k, unsigned numHashes)
		: m_seq(seq), m_seqLen(strlen(seq)), m_k(k),
		m_numHashes(numHashes), m_pos(0), m_hash1(0), m_rcHash1(0)
	{
		assert(seq != NULL);
		m_nextInvalidChar = strspn(seq, ACGT_CHARS);
		next();
	}

	/** get reference to hash values for current k-mer */
	const std::vector<hash_t>& operator*() const
	{
		assert(m_pos + m_k <= m_seqLen);
		return m_hashes;
	}

	/** get pointer to hash values for current k-mer */
	const std::vector<hash_t>* operator->() const
	{
		assert(m_pos + m_k <= m_seqLen);
		return &m_hashes;
	}

	/** test equality with another iterator */
	bool operator==(const RollingHashIterator& it) const
	{
		return m_pos == it.m_pos;
	}

	/** test inequality with another iterator */
	bool operator!=(const RollingHashIterator& it) const
	{
		return !(*this == it);
	}

	/** pre-increment operator */
	RollingHashIterator& operator++()
	{
		assert(m_pos + m_k <= m_seqLen);
		++m_pos;
		next();
		return *this;
	}

	/** post-increment operator */
	RollingHashIterator operator++(int)
	{
		RollingHashIterator it = *this;
		++*this;
		return it;
	}

    /** iterator pointing to one past last element */
	static const RollingHashIterator end()
	{
		return RollingHashIterator();
	}

	/** return position of current k-mer */
	unsigned pos() const
	{
		return m_pos;
	}

private:

	/** DNA sequence being hashed */
	const char* m_seq;
	/** length of DNA sequence being hashed */
	const unsigned m_seqLen;
	/** k-mer size */
	unsigned m_k;
	/** number of hash values to compute for each k-mer */
	unsigned m_numHashes;
	/** position of current k-mer */
	size_t m_pos;
	/** position of next non-ACGT char */
	size_t m_nextInvalidChar;
	/** hash values for the current k-mer */
	std::vector<hash_t> m_hashes;
	/** Value of first hash function on current k-mer. */
	hash_t m_hash1;
	/**
	 * Value of first hash function on reverse complement
	 * of current k-mer */
	hash_t m_rcHash1;
};

#endif
