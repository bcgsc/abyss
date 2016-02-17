#ifndef MASKED_KMER_H
#define MASKED_KMER_H 1

#include "Common/Kmer.h"
#include "Common/Hash.h"
#include "Common/Sequence.h"
#include <iostream>

class MaskedKmer : public Kmer
{
public:
	friend class Kmer;

	/** Default constructor */
	MaskedKmer() : Kmer() {}

	/**
	 * Constructor.
	 * @param seq k-mer sequence
	 */
	explicit MaskedKmer(const Sequence& seq) : Kmer(seq)
	{
		std::string kmerMask(length(), '1');
		initKmerMask(kmerMask);
	}

	/**
	 * Constructor.
	 * @param seq k-mer sequence
	 * @param kmerMask bitmask of "don't care" positions
	 */
	MaskedKmer(const Sequence& seq, const std::string& kmerMask) : Kmer(seq)
	{
		initKmerMask(kmerMask);
	}

	/** Compare this k-mer to another */
	int compare(const Kmer& other) const
	{
		if (m_nonTrivialMask) {
			Kmer kmer1(*this), kmer2(other);
			maskKmer(kmer1);
			maskKmer(kmer2);
			return kmer1.compare(kmer2);
		}
		return Kmer::compare(other);
	}

	/** Equality operator */
	bool operator==(const MaskedKmer& other) const
	{
		return compare(other) == 0;
	}

	/** Inequality operator */
	bool operator!=(const MaskedKmer& other) const
	{
		return compare(other) != 0;
	}

	/** Less-than operator */
	bool operator<(const MaskedKmer& other) const
	{
		return compare(other) < 0;
	}

protected:

	/** Initialize variables related to masking of "don't care" positions */
	void initKmerMask(const std::string& kmerMask)
	{
		/* check that mask has length k */
		assert(kmerMask.length() == length());
		m_kmerMask = kmerMask;

		/* use more efficient unmasked operations if k-mer mask is all "1"s */
		m_nonTrivialMask = false;
		for (size_t i = 0; i < m_kmerMask.length(); ++i) {
			if (m_kmerMask.at(i) == '0') {
				m_nonTrivialMask = true;
				break;
			}
		}
	}

	/** Mask out don't care positions by changing them to 'A' */
	void maskKmer(Kmer& kmer) const
	{
		assert(m_kmerMask.length() == length());
		for(size_t i = 0; i < m_kmerMask.length(); ++i) {
			if (m_kmerMask.at(i) == '0')
				kmer.set(i, baseToCode('A'));
		}
	}

	/* true if kmerMask contains one or more '0's */
	bool m_nonTrivialMask;
	/* bitmask indicating "don't care" positions */
	std::string m_kmerMask;
};

#endif
