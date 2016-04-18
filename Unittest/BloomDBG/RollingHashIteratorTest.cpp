#include "BloomDBG/RollingHashIterator.h"

#include <gtest/gtest.h>
#include <string>

using namespace std;

TEST(RollingHashIterator, reverseComplement)
{
	const unsigned k = 6;
	const unsigned numHashes = 1;
	const char* seq = "GCAATGT";
	const char* rcSeq = "ACATTGC";

	/** hash forward sequence */

	RollingHashIterator it(seq, numHashes, k);
	size_t kmer1Hash, kmer2Hash;
	kmer1Hash = (*it)[0];
	++it;
	kmer2Hash = (*it)[0];
	++it;
	ASSERT_EQ(RollingHashIterator::end(), it);

	/** hash reverse complement sequence */

	RollingHashIterator rcIt(rcSeq, numHashes, k);
	size_t rcKmer1Hash, rcKmer2Hash;
	rcKmer2Hash = (*rcIt)[0];
	++rcIt;
	rcKmer1Hash = (*rcIt)[0];
	++rcIt;
	ASSERT_EQ(RollingHashIterator::end(), rcIt);

	/** check hash values are the same for forward and reverse complement */

	ASSERT_EQ(kmer1Hash, rcKmer1Hash);
	ASSERT_EQ(kmer2Hash, rcKmer2Hash);
}

TEST(RollingHashIterator, badKmers)
{
	const unsigned k = 3;
	const unsigned numHashes = 1;

    /* skip bad k-mers in middle of sequence */

	const char* seq = "AAANAAA";
    RollingHashIterator it(seq, numHashes, k);
	ASSERT_EQ(0u, it.pos());
	++it;
	ASSERT_EQ(4u, it.pos());
	++it;
	ASSERT_EQ(RollingHashIterator::end(), it);

	/* all bad k-mers */

	const char* seq2 = "NNNNNNN";
	RollingHashIterator it2(seq2, numHashes, k);
	ASSERT_EQ(RollingHashIterator::end(), it2);
}

TEST(RollingHashIterator, seqShorterThanK)
{
	const unsigned k = 5;
	const unsigned numHashes = 1;
	const char* seq = "ACGT";

	RollingHashIterator it(seq, numHashes, k);
	ASSERT_EQ(RollingHashIterator::end(), it);
}

TEST(RollingHashIterator, emptySeq)
{
	const unsigned k = 3;
	const unsigned numHashes = 1;
	const char* seq = "";

	RollingHashIterator it(seq, numHashes, k);
	ASSERT_EQ(RollingHashIterator::end(), it);
}

TEST(RollingHashIterator, spacedSeed)
{
	const unsigned k = 5;
	const unsigned numHashes = 1;
	const char* seq = "AGNNGC";
	const char* rcSeq = "GCNNCT";
	Kmer::setLength(k);
	MaskedKmer::setMask("10001");

	/** hash forward sequence */

	RollingHashIterator it(seq, numHashes, k);
	size_t kmer1Hash, kmer2Hash;
	kmer1Hash = (*it)[0];
	++it;
	kmer2Hash = (*it)[0];
	++it;
	ASSERT_EQ(RollingHashIterator::end(), it);

	/** hash reverse complement sequence */

	RollingHashIterator rcIt(rcSeq, numHashes, k);
	size_t rcKmer1Hash, rcKmer2Hash;
	rcKmer2Hash = (*rcIt)[0];
	++rcIt;
	rcKmer1Hash = (*rcIt)[0];
	++rcIt;
	ASSERT_EQ(RollingHashIterator::end(), rcIt);

	/** check hash values are the same for forward and reverse complement */

	ASSERT_EQ(kmer1Hash, rcKmer1Hash);
	ASSERT_EQ(kmer2Hash, rcKmer2Hash);
}
