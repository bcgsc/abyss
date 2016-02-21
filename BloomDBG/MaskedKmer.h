#ifndef MASKED_KMER_H
#define MASKED_KMER_H 1

#include "Common/Kmer.h"
#include "Common/Hash.h"
#include "Common/Sequence.h"
#include <iostream>
#include <string>
#include <cstdlib>

class MaskedKmer : public Kmer
{
public:

	/** Default constructor */
	MaskedKmer() : Kmer() {}

	/**
	 * Constructor.
	 * @param seq k-mer sequence
	 */
	explicit MaskedKmer(const Sequence& seq) : Kmer(seq) {}

	/** Set global k-mer mask (a.k.a. spaced seed) */
	static void setMask(const std::string& kmerMask)
	{
		/* setLength() must be called before setMask() */
		assert(length() > 0);

		/* set global bitmask */
		mask() = kmerMask;

		/* empty mask is equivalent to string of '1's */
		if (kmerMask.empty())
			return;

		/* check for valid spaced seed pattern */
		if (mask().length() != length()) {
			std::cerr << "error: spaced seed must be exactly k bits long\n";
			exit(EXIT_FAILURE);
		} else if (mask().find_first_not_of("01") != std::string::npos) {
			std::cerr << "error: spaced seed must contain only '0's or '1's\n";
			exit(EXIT_FAILURE);
		} else if (*mask().begin() != '1' || *mask().rbegin() != '1') {
			std::cerr << "error: spaced seed must begin and end with '1's\n";
			exit(EXIT_FAILURE);
		}
	}

	/** Get global k-mer mask */
	static std::string& mask()
	{
		static std::string s_kmerMask;
		return s_kmerMask;
	}

	/** Compare this k-mer to another */
	int compare(const Kmer& other) const
	{
		if (!mask().empty()) {
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

	/** Mask out don't care positions by changing them to 'A' */
	static void maskKmer(Kmer& kmer)
	{
		if (mask().empty())
			return;

		assert(mask().length() == length());
		for(size_t i = 0; i < mask().length(); ++i) {
			if (mask().at(i) == '0')
				kmer.set(i, baseToCode('A'));
		}
	}
};

/** Return the reverse complement of the specified k-mer. */
static inline MaskedKmer reverseComplement(const MaskedKmer& seq)
{
	MaskedKmer rc(seq);
	rc.reverseComplement();
	return rc;
}

/** Define default hash function for use with STL containers */
NAMESPACE_STD_HASH_BEGIN
template <> struct hash<MaskedKmer> {
	size_t operator()(const MaskedKmer& kmer) const
	{
		MaskedKmer kmerCopy(kmer);
		MaskedKmer::maskKmer(kmerCopy);
		return kmerCopy.getHashCode();
	}
};
NAMESPACE_STD_HASH_END

#endif
