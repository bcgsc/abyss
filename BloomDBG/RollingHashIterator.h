#ifndef ROLLING_HASH_ITERATOR_H
#define ROLLING_HASH_ITERATOR_H 1

#include <cstring>
#include <vector>
#include <cassert>
#include <limits>
#include <string>
#include <algorithm>
#include <cctype>
#include <deque>
#include "BloomDBG/RollingHash.h"

/**
 * Permitted characters in k-mers. All k-mers containing
 * other characters will be skipped.
 */
#define ACGT_CHARS "ACGT"

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
private:

	/**
	 * Advance iterator right to the next valid k-mer.
	 */
	void next()
	{
		if (m_seq.length() < m_k) {
			m_pos = std::numeric_limits<std::size_t>::max();
			return;
		}

		const std::string& spacedSeed = MaskedKmer::mask();

		while(m_pos < m_seq.length() - m_k + 1) {

			/* skip k-mers with non-ACGT chars in unmasked positions */

			while (!m_badCharPos.empty() && m_badCharPos.front() < m_pos)
				m_badCharPos.pop_front();

			if (!m_badCharPos.empty() && m_badCharPos.front() < m_pos + m_k) {
				/* empty spaced seed is equivalent to a string of '1's */
				if (spacedSeed.empty()) {
					m_rollNextHash = false;
					m_pos = m_badCharPos.front() + 1;
					continue;
				}
				bool goodKmer = true;
				assert(spacedSeed.length() == m_k);
				for (size_t i = 0; i < m_badCharPos.size() &&
					m_badCharPos.at(i) < m_pos + m_k; ++i) {
					size_t kmerPos = m_badCharPos.at(i) - m_pos;
					if (spacedSeed.at(kmerPos) == '1') {
						goodKmer = false;
						break;
					}
				}
				if (!goodKmer) {
					m_rollNextHash = false;
					++m_pos;
					continue;
				}
			}

			/* we are positioned at the next valid k-mer */

			if (!m_rollNextHash) {
				/* we don't have hash values for the
				 * preceding k-mer, so we must compute
				 * the hash values from scratch */
				m_rollingHash.reset(m_seq.substr(m_pos, m_k));
				m_rollNextHash = true;
			} else {
				/* compute new hash values based on
				 * hash values of preceding k-mer */
				assert(m_pos > 0);
				m_rollingHash.rollRight(m_seq.c_str() + m_pos - 1,
					m_seq[m_pos + m_k - 1]);
			}
			m_rollingHash.getHashes(m_hashes);
			return;

		}

		/* there are no more valid k-mers */
		m_pos = std::numeric_limits<std::size_t>::max();
	}

	/**
	 * Advance iterator left to the prev valid k-mer.
	 */
	void prev()
	{
		if (m_seq.length() < m_k) {
			m_pos = std::numeric_limits<std::size_t>::max();
			return;
		}

		const std::string& spacedSeed = MaskedKmer::mask();

		//while(m_pos < m_seq.length() - m_k + 1) {
		while(m_pos < m_seq.length() - m_k + 1) {

			/* skip k-mers with non-ACGT chars in unmasked positions */

			//TODO
			//while (!m_badCharPos.empty() && m_badCharPos.front() < m_pos)
			//	m_badCharPos.pop_front();
			while (!m_badCharPos.empty() && m_badCharPos.back() > m_pos + m_k - 1)
				m_badCharPos.pop_back();

			//if (!m_badCharPos.empty() && m_badCharPos.front() < m_pos + m_k) {
			if (!m_badCharPos.empty() && m_badCharPos.back() >= m_pos) {
				/* empty spaced seed is equivalent to a string of '1's */
				if (spacedSeed.empty()) {
					m_rollNextHash = false;
					//m_pos = m_badCharPos.front() + 1;
					m_pos = m_badCharPos.back() - m_k;
					continue;
				}
				bool goodKmer = true;
				assert(spacedSeed.length() == m_k);
				//for (size_t i = 0; i < m_badCharPos.size() &&
				//	m_badCharPos.at(i) < m_pos + m_k; ++i) {
				//	size_t kmerPos = m_badCharPos.at(i) - m_pos;
				//	if (spacedSeed.at(kmerPos) == '1') {
				//		goodKmer = false;
				//		break;
				//	}
				//}
				for (size_t i = m_badCharPos.size() - 1; i != std::numeric_limits<std::size_t>::max() &&
				    m_badCharPos.at(i) >= m_pos; --i) {
					size_t kmerPos = m_badCharPos.at(i) - m_pos;
					if (spacedSeed.at(kmerPos) == '1') {
						goodKmer = false;
						break;
					}
				}
				if (!goodKmer) {
					m_rollNextHash = false;
					--m_pos;
					continue;
				}
			}

			/* we are positioned at the next valid k-mer */

			if (!m_rollNextHash) {
				/* we don't have hash values for the
				 * preceding k-mer, so we must compute
				 * the hash values from scratch */
				//std::cerr << m_pos << std::endl;
				m_rollingHash.reset(m_seq.substr(m_pos, m_k));
				m_rollNextHash = true;
			} else {
				/* compute new hash values based on
				 * hash values of preceding k-mer */
				assert(m_pos < m_seq.length() - m_k);
				m_rollingHash.rollLeft(m_seq[m_pos], m_seq.c_str() + m_pos + 1);
			}
			m_rollingHash.getHashes(m_hashes);
			return;

		}

		/* there are no more valid k-mers */
		m_pos = std::numeric_limits<std::size_t>::max();
	}

public:

	typedef uint64_t hash_t;
	/**
	 * Default constructor. Creates an iterator pointing to
	 * the end of the iterator range.
	 */
	RollingHashIterator() : m_numHashes(0), m_k(0),
		m_rollingHash(m_numHashes, m_k),
		m_pos(std::numeric_limits<std::size_t>::max()) {}

	/**
	 * Constructor.
	 * @param seq DNA sequence to be hashed
	 * @param k k-mer size
	 * for each k-mer
	 */
	RollingHashIterator(const std::string& seq, unsigned numHashes, unsigned k, bool back = false)
		: m_seq(seq), m_numHashes(numHashes), m_k(k),
		m_rollingHash(m_numHashes, m_k), m_rollNextHash(false)
	{
		if (!back) {
			m_pos = 0;
			init();
		} else {
			m_pos = m_seq.length() - m_k;
			initFromBack();
		}
	}

	/**
	 * Initialize internal state of iterator.
	 */
	void init()
	{
		/* note: empty spaced seed indicates no masking (string of '1's) */
		assert(MaskedKmer::mask().empty() || MaskedKmer::mask().length() == m_k);

		/* convert sequence to upper case */
		std::transform(m_seq.begin(), m_seq.end(), m_seq.begin(), ::toupper);

		/* record positions of non-ACGT chars */
		size_t i = m_seq.find_first_not_of(ACGT_CHARS);
		while (i != std::string::npos) {
			m_badCharPos.push_back(i);
			i = m_seq.find_first_not_of(ACGT_CHARS, i + 1);
		}

		/* find first "good" k-mer in sequence */
		next();
	}

	/**
	 * Initialize internal state of iterator.
	 */
	void initFromBack()
	{
		/* note: empty spaced seed indicates no masking (string of '1's) */
		assert(MaskedKmer::mask().empty() || MaskedKmer::mask().length() == m_k);

		/* convert sequence to upper case */
		std::transform(m_seq.begin(), m_seq.end(), m_seq.begin(), ::toupper);

		/* record positions of non-ACGT chars */
		size_t i = m_seq.find_first_not_of(ACGT_CHARS);
		while (i != std::string::npos) {
			m_badCharPos.push_back(i);
			i = m_seq.find_first_not_of(ACGT_CHARS, i + 1);
		}
		//std::cerr << "num bad bases" <<std::endl;
		//std::cerr << m_badCharPos.size() <<std::endl;

		/* find first "good" k-mer in sequence */
		prev();
	}

	/** get reference to hash values for current k-mer */
	const hash_t* operator*() const
	{
		assert(m_pos + m_k <= m_seq.length());
		return m_hashes;
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
		++m_pos;
		next();
		return *this;
	}
	
	/** pre-decrement operator */
	RollingHashIterator& operator--()
	{
		//std::cerr << m_pos << std::endl;
		--m_pos;
		prev();
		return *this;
	}

	/** post-increment operator */
	RollingHashIterator operator++(int)
	{
		RollingHashIterator it = *this;
		++*this;
		return it;
	}

	/** post-decrement operator */
	RollingHashIterator operator--(int)
	{
		RollingHashIterator it = *this;
		--*this;
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

	/** return k-mer at current position */
	std::string kmer(bool mask=false) const
	{
		std::string kmer(m_seq, m_pos, m_k);
		const std::string& spacedSeed = MaskedKmer::mask();
		if (mask && !spacedSeed.empty()) {
			assert(spacedSeed.length() == m_k);
			for(size_t i = 0; i < spacedSeed.length(); ++i) {
				if (spacedSeed.at(i) == '0')
					kmer.at(i) = 'N';
			}
		}
		return kmer;
	}

	/** return RollingHash object for current state */
	RollingHash rollingHash()
	{
		return m_rollingHash;
	}

private:

	/** DNA sequence being hashed */
	std::string m_seq;
	/** number of hash functions */
	unsigned m_numHashes;
	/** hash values */
	hash_t m_hashes[MAX_HASHES];
	/** k-mer size */
	unsigned m_k;
	/** internal state for rolling hash */
	RollingHash m_rollingHash;
	/** true whenever we can "roll" the hash values for
	 * the current k-mer to compute the hash values for the
	 * next k-mer */
	bool m_rollNextHash;
	/** position of current k-mer */
	size_t m_pos;
	/** positions of non-ACGT chars in sequence */
	std::deque<size_t> m_badCharPos;
};

#endif
