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
	const dynamic_bitset<> m_kmerMask;

	RollingHashTest() : m_numHashes(2), m_k(4), m_kmerMask(string("1101")) {}
};

TEST_F(RollingHashTest, kmerMask)
{
	RollingHash kmer1Hash("GACG", m_numHashes, m_k, m_kmerMask);
	RollingHash kmer2Hash("GATG", m_numHashes, m_k, m_kmerMask);
	ASSERT_EQ(kmer1Hash, kmer2Hash);
}

TEST_F(RollingHashTest, rollRight)
{
	RollingHash middleKmerHash("ACGT", m_numHashes, m_k);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k);

	middleKmerHash.rollRight('A', 'C');
	ASSERT_EQ(rightKmerHash.getHash(), middleKmerHash.getHash());
}

TEST_F(RollingHashTest, rollRightMasked)
{
	RollingHash middleKmerHash("ACGG", m_numHashes, m_k, m_kmerMask);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k, m_kmerMask);

	middleKmerHash.rollRight('A', 'C');
	ASSERT_EQ(rightKmerHash.getHash(), middleKmerHash.getHash());
}

TEST_F(RollingHashTest, rollLeft)
{
	RollingHash leftKmerHash("GACG", m_numHashes, m_k);
	RollingHash middleKmerHash("ACGT", m_numHashes, m_k);

	middleKmerHash.rollLeft('G', 'T');
	ASSERT_EQ(leftKmerHash.getHash(), middleKmerHash.getHash());
}

TEST_F(RollingHashTest, rollLeftMasked)
{
	RollingHash leftKmerHash("GACG", m_numHashes, m_k, m_kmerMask);
	RollingHash middleKmerHash("AGGT", m_numHashes, m_k, m_kmerMask);

	middleKmerHash.rollLeft('G', 'T');
	ASSERT_EQ(leftKmerHash.getHash(), middleKmerHash.getHash());
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
	RollingHash middleKmerHash("ACGG", m_numHashes, m_k, m_kmerMask);
	RollingHash rightKmerHash("CGTC", m_numHashes, m_k, m_kmerMask);

	middleKmerHash.reset("CGGC");
	ASSERT_EQ(rightKmerHash.getHash(), middleKmerHash.getHash());
}
