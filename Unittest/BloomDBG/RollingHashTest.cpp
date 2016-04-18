#include "BloomDBG/RollingHash.h"

#include <gtest/gtest.h>
#include <string>
#include <algorithm>
#include "boost/dynamic_bitset.hpp"

using namespace std;
using namespace boost;

/** test fixture for RollingHash tests */
class RollingHashTest : public ::testing::Test
{
protected:

	const unsigned m_numHashes;
	const unsigned m_k;
	const string m_kmerMask;

	RollingHashTest() : m_numHashes(2), m_k(4)
	{
		Kmer::setLength(m_k);
	}
};

TEST_F(RollingHashTest, kmerMask)
{
	MaskedKmer::setMask("1001");
	RollingHash kmer1Hash("GCCG", m_numHashes, m_k);
	RollingHash kmer2Hash("GTTG", m_numHashes, m_k);
	ASSERT_EQ(kmer1Hash, kmer2Hash);
}

TEST_F(RollingHashTest, rollRight)
{
	MaskedKmer::mask().clear();
	RollingHash leftKmerHash("GACG", m_numHashes, m_k);
	RollingHash middleKmerHash("ACGT", m_numHashes, m_k);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k);

	leftKmerHash.rollRight("GACG", 'T');
	ASSERT_EQ(middleKmerHash, leftKmerHash);
	leftKmerHash.rollRight("ACGT", 'C');
	ASSERT_EQ(rightKmerHash, leftKmerHash);
}

TEST_F(RollingHashTest, rollRightMasked)
{
	MaskedKmer::setMask("1001");
	RollingHash leftKmerHash("GACG", m_numHashes, m_k);
	RollingHash middleKmerHash("ACGT", m_numHashes, m_k);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k);

	leftKmerHash.rollRight("GACG", 'T');
	ASSERT_EQ(middleKmerHash, leftKmerHash);
	leftKmerHash.rollRight("ACGT", 'C');
	ASSERT_EQ(rightKmerHash, leftKmerHash);
}

TEST_F(RollingHashTest, rollRightMaskedMismatch)
{
	MaskedKmer::setMask("1001");

	const char* origSeq    = "GACGTC";
	const char* mutatedSeq = "GACTTC";

	RollingHash left(origSeq, m_numHashes, m_k);
	RollingHash middle(origSeq + 1, m_numHashes, m_k);
	RollingHash right(origSeq + 2, m_numHashes, m_k);

	RollingHash mutated(mutatedSeq, m_numHashes, m_k);

	ASSERT_NE(left, mutated);
	mutated.rollRight(mutatedSeq, 'T');
	ASSERT_EQ(middle, mutated);
	mutated.rollRight(mutatedSeq + 1, 'C');
	ASSERT_EQ(right, mutated);
}

TEST_F(RollingHashTest, rollLeft)
{
	MaskedKmer::mask().clear();

	RollingHash leftKmerHash("GACG", m_numHashes, m_k);
	RollingHash middleKmerHash("ACGT", m_numHashes, m_k);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k);

	rightKmerHash.rollLeft('A', "CGTC");
	ASSERT_EQ(middleKmerHash, rightKmerHash);
	rightKmerHash.rollLeft('G', "ACGT");
	ASSERT_EQ(leftKmerHash, rightKmerHash);
}

TEST_F(RollingHashTest, rollLeftMasked)
{
	MaskedKmer::setMask("1001");

	RollingHash leftKmerHash("GACG", m_numHashes, m_k);
	RollingHash middleKmerHash("ACGT", m_numHashes, m_k);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k);

	rightKmerHash.rollLeft('A', "CGTC");
	ASSERT_EQ(middleKmerHash, rightKmerHash);
	rightKmerHash.rollLeft('G', "ACGT");
	ASSERT_EQ(leftKmerHash, rightKmerHash);
}

TEST_F(RollingHashTest, rollLeftMaskedMismatch)
{
	MaskedKmer::setMask("1001");

	const char* origSeq    = "GACGTC";
	const char* mutatedSeq = "GAGGTC";

	RollingHash left(origSeq, m_numHashes, m_k);
	RollingHash middle(origSeq + 1, m_numHashes, m_k);
	RollingHash right(origSeq + 2, m_numHashes, m_k);

	RollingHash mutated(mutatedSeq + 2, m_numHashes, m_k);

	ASSERT_NE(right, mutated);
	mutated.rollLeft('A', mutatedSeq + 2);
	ASSERT_EQ(middle, mutated);
	mutated.rollLeft('G', mutatedSeq + 1);
	ASSERT_EQ(left, mutated);
}

TEST_F(RollingHashTest, reset)
{
	MaskedKmer::mask().clear();

	RollingHash middleKmerHash("ACGT", m_numHashes, m_k);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k);

	middleKmerHash.reset("CGTC");
	ASSERT_EQ(rightKmerHash, middleKmerHash);
}

TEST_F(RollingHashTest, resetMasked)
{
	MaskedKmer::setMask("1001");

	RollingHash middleKmerHash("ACGT", m_numHashes, m_k);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k);

	/*
	 * Note: third base of middleKmerHash is intentionally set to 'G'
	 * instead of 'T'. However, the hash values should
	 * still match the rightKmerHash due to the effect of
	 * the k-mer mask.
	 */
	middleKmerHash.reset("CGGC");
	ASSERT_EQ(rightKmerHash, middleKmerHash);
}

TEST_F(RollingHashTest, setBase)
{
	MaskedKmer::mask().clear();

	char kmer1[] = "ACGT";
	char kmer2[] = "ACCT";

	RollingHash hash1(kmer1, m_numHashes, m_k);
	RollingHash hash2(kmer2, m_numHashes, m_k);

	ASSERT_NE(hash2, hash1);
	hash1.setBase(kmer1, 2, 'C');
	ASSERT_EQ(0, strcmp(kmer1, kmer2));
	ASSERT_EQ(hash2, hash1);
}

TEST_F(RollingHashTest, setBaseMasked)
{
	MaskedKmer::setMask("1101");

	char kmer1[] = "ACGT";
	char kmer2[] = "ACCT";

	RollingHash hash1(kmer1, m_numHashes, m_k);
	RollingHash hash2(kmer2, m_numHashes, m_k);

	/* hashes should agree since mismatch is in masked position */
	ASSERT_EQ(hash2, hash1);
	ASSERT_NE(0, strcmp(kmer1, kmer2));

	/* fix mismatch in masked position (hash values shouldn't change) */
	hash1.setBase(kmer1, 2, 'C');
	ASSERT_EQ(hash2, hash1);
	ASSERT_EQ(0, strcmp(kmer1, kmer2));

	/* create mismatch in unmasked position (hash value should now differ) */
	hash1.setBase(kmer1, 1, 'G');
	ASSERT_NE(hash2, hash1);
	ASSERT_NE(0, strcmp(kmer1, kmer2));
}
