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

	const unsigned m_k;
	const string m_kmerMask;

	RollingHashTest() : m_k(4)
	{
		Kmer::setLength(m_k);
		MaskedKmer::setMask("1001");
	}
};

TEST_F(RollingHashTest, kmerMask)
{
	RollingHash kmer1Hash("GCCG", m_k);
	RollingHash kmer2Hash("GTTG", m_k);
	ASSERT_EQ(kmer1Hash, kmer2Hash);
}

TEST_F(RollingHashTest, rollRight)
{
	RollingHash leftKmerHash("GACG", m_k);
	RollingHash middleKmerHash("ACGT", m_k);
	RollingHash rightKmerHash("CGTC", m_k);

	leftKmerHash.rollRight("GACG", "ACGT");
	ASSERT_EQ(middleKmerHash, leftKmerHash);
	leftKmerHash.rollRight("ACGT", "CGTC");
	ASSERT_EQ(rightKmerHash, leftKmerHash);
}

TEST_F(RollingHashTest, rollRightMasked)
{
	RollingHash leftKmerHash("GACG", m_k);
	RollingHash middleKmerHash("ACGT", m_k);
	RollingHash rightKmerHash("CGTC", m_k);

	leftKmerHash.rollRight("GACG", "ACGT");
	ASSERT_EQ(middleKmerHash, leftKmerHash);
	leftKmerHash.rollRight("ACGT", "CGTC");
	ASSERT_EQ(rightKmerHash, leftKmerHash);
}

TEST_F(RollingHashTest, rollRightMaskedMismatch)
{
	/*
	 * Note: last base of leftKmerHash is intentionally set to 'T'
	 * instead of 'G'. However, the rolled hash values should
	 * still match the middle/right k-mers due to the effect of
	 * the k-mer mask.
	 */
	RollingHash leftKmerHash("GACT", m_k);
	RollingHash middleKmerHash("ACGT", m_k);
	RollingHash rightKmerHash("CGTC", m_k);

	leftKmerHash.rollRight("GACT", "ACGT");
	ASSERT_EQ(middleKmerHash, leftKmerHash);
	leftKmerHash.rollRight("ACGT", "CGTC");
	ASSERT_EQ(rightKmerHash, leftKmerHash);
}

TEST_F(RollingHashTest, rollLeft)
{
	RollingHash leftKmerHash("GACG", m_k);
	RollingHash middleKmerHash("ACGT", m_k);
	RollingHash rightKmerHash("CGTC", m_k);

	rightKmerHash.rollLeft("ACGT", "CGTC");
	ASSERT_EQ(middleKmerHash, rightKmerHash);
	rightKmerHash.rollLeft("GACG", "ACGT");
	ASSERT_EQ(leftKmerHash, rightKmerHash);
}

TEST_F(RollingHashTest, rollLeftMasked)
{
	RollingHash leftKmerHash("GACG", m_k);
	RollingHash middleKmerHash("ACGT", m_k);
	RollingHash rightKmerHash("CGTC", m_k);

	rightKmerHash.rollLeft("ACGT", "CGTC");
	ASSERT_EQ(middleKmerHash, rightKmerHash);
	rightKmerHash.rollLeft("GACG", "ACGT");
	ASSERT_EQ(leftKmerHash, rightKmerHash);
}

TEST_F(RollingHashTest, rollLeftMaskedMismatch)
{
	RollingHash leftKmerHash("GACG", m_k);
	RollingHash middleKmerHash("ACGT", m_k);
	/*
	 * Note: first base of rightKmerHash is intentionally set to 'G'
	 * instead of 'C'. However, the rolled hash values should
	 * still match the left/middle k-mers due to the effect of
	 * the k-mer mask.
	 */
	RollingHash rightKmerHash("GGTC", m_k);

	rightKmerHash.rollLeft("ACGT", "GGTC");
	ASSERT_EQ(middleKmerHash, rightKmerHash);
	rightKmerHash.rollLeft("GACG", "ACGT");
	ASSERT_EQ(leftKmerHash, rightKmerHash);
}

TEST_F(RollingHashTest, reset)
{
	RollingHash middleKmerHash("ACGT", m_k);
	RollingHash rightKmerHash("CGTC", m_k);

	middleKmerHash.reset("CGTC");
	ASSERT_EQ(rightKmerHash, middleKmerHash);
}

TEST_F(RollingHashTest, resetMasked)
{
	RollingHash middleKmerHash("ACGT", m_k);
	RollingHash rightKmerHash("CGTC", m_k);

	/*
	 * Note: third base of middleKmerHash is intentionally set to 'G'
	 * instead of 'T'. However, the hash values should
	 * still match the rightKmerHash due to the effect of
	 * the k-mer mask.
	 */
	middleKmerHash.reset("CGGC");
	ASSERT_EQ(rightKmerHash, middleKmerHash);
}
