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

	RollingHashTest() : m_numHashes(2), m_k(4), m_kmerMask("1001") {}
};

TEST_F(RollingHashTest, kmerMask)
{
	RollingHash kmer1Hash("GCCG", m_numHashes, m_k, m_kmerMask);
	RollingHash kmer2Hash("GTTG", m_numHashes, m_k, m_kmerMask);
	ASSERT_EQ(kmer1Hash, kmer2Hash);
}

TEST_F(RollingHashTest, rollRight)
{
	RollingHash leftKmerHash("GACG", m_numHashes, m_k);
	RollingHash middleKmerHash("ACGT", m_numHashes, m_k);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k);

	leftKmerHash.rollRight('G', 'T');
	ASSERT_EQ(middleKmerHash.getHash(), leftKmerHash.getHash());
	leftKmerHash.rollRight('A', 'C');
	ASSERT_EQ(rightKmerHash.getHash(), leftKmerHash.getHash());
}

TEST_F(RollingHashTest, rollRightMasked)
{
	RollingHash leftKmerHash("GACG", m_numHashes, m_k, m_kmerMask);
	RollingHash middleKmerHash("ACGT", m_numHashes, m_k, m_kmerMask);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k, m_kmerMask);

	leftKmerHash.rollRight('G', 'T');
	ASSERT_EQ(middleKmerHash.getHash(), leftKmerHash.getHash());
	leftKmerHash.rollRight('A', 'C');
	ASSERT_EQ(rightKmerHash.getHash(), leftKmerHash.getHash());
}

TEST_F(RollingHashTest, rollRightMaskedMismatch)
{
	/*
	 * Note: last base of leftKmerHash is intentionally set to 'T'
	 * instead of 'G'. However, the rolled hash values should
	 * still match the middle/right k-mers due to the effect of
	 * the k-mer mask.
	 */
	RollingHash leftKmerHash("GACT", m_numHashes, m_k, m_kmerMask);
	RollingHash middleKmerHash("ACGT", m_numHashes, m_k, m_kmerMask);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k, m_kmerMask);

	leftKmerHash.rollRight('G', 'T');
	ASSERT_EQ(middleKmerHash.getHash(), leftKmerHash.getHash());
	leftKmerHash.rollRight('A', 'C');
	ASSERT_EQ(rightKmerHash.getHash(), leftKmerHash.getHash());
}

TEST_F(RollingHashTest, rollLeft)
{
	RollingHash leftKmerHash("GACG", m_numHashes, m_k);
	RollingHash middleKmerHash("ACGT", m_numHashes, m_k);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k);

	rightKmerHash.rollLeft('A', 'C');
	ASSERT_EQ(middleKmerHash.getHash(), rightKmerHash.getHash());
	rightKmerHash.rollLeft('G', 'T');
	ASSERT_EQ(leftKmerHash.getHash(), rightKmerHash.getHash());
}

TEST_F(RollingHashTest, rollLeftMasked)
{
	RollingHash leftKmerHash("GACG", m_numHashes, m_k, m_kmerMask);
	RollingHash middleKmerHash("ACGT", m_numHashes, m_k, m_kmerMask);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k, m_kmerMask);

	rightKmerHash.rollLeft('A', 'C');
	ASSERT_EQ(middleKmerHash.getHash(), rightKmerHash.getHash());
	rightKmerHash.rollLeft('G', 'T');
	ASSERT_EQ(leftKmerHash.getHash(), rightKmerHash.getHash());
}

TEST_F(RollingHashTest, rollLeftMaskedMismatch)
{
	RollingHash leftKmerHash("GACG", m_numHashes, m_k, m_kmerMask);
	RollingHash middleKmerHash("ACGT", m_numHashes, m_k, m_kmerMask);
	/*
	 * Note: first base of rightKmerHash is intentionally set to 'G'
	 * instead of 'C'. However, the rolled hash values should
	 * still match the left/middle k-mers due to the effect of
	 * the k-mer mask.
	 */
	RollingHash rightKmerHash("GGTC", m_numHashes, m_k, m_kmerMask);

	rightKmerHash.rollLeft('A', 'C');
	ASSERT_EQ(middleKmerHash.getHash(), rightKmerHash.getHash());
	rightKmerHash.rollLeft('G', 'T');
	ASSERT_EQ(leftKmerHash.getHash(), rightKmerHash.getHash());
}

TEST_F(RollingHashTest, reset)
{
	RollingHash middleKmerHash("ACGT", m_numHashes, m_k);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k);

	middleKmerHash.reset("CGTC");
	ASSERT_EQ(rightKmerHash.getHash(), middleKmerHash.getHash());
}

TEST_F(RollingHashTest, resetMasked)
{
	RollingHash middleKmerHash("ACGT", m_numHashes, m_k, m_kmerMask);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k, m_kmerMask);

	/*
	 * Note: third base of middleKmerHash is intentionally set to 'G'
	 * instead of 'T'. However, the hash values should
	 * still match the rightKmerHash due to the effect of
	 * the k-mer mask.
	 */
	middleKmerHash.reset("CGGC");
	ASSERT_EQ(rightKmerHash.getHash(), middleKmerHash.getHash());
}
