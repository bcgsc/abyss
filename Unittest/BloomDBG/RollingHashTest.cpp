#include "BloomDBG/RollingHash.h"

#include <gtest/gtest.h>

using namespace std;

/** test fixture for RollingHash tests */
class RollingHashTest : public ::testing::Test
{
protected:

	const unsigned m_numHashes;
	const unsigned m_k;
	const string m_leftKmer;
	const string m_middleKmer;
	const string m_rightKmer;

	RollingHashTest() : m_numHashes(2), m_k(4),
		m_leftKmer("GACG"), m_middleKmer("ACGT"), m_rightKmer("CGTC") {}
};

TEST_F(RollingHashTest, peekRight)
{
	RollingHash middleKmerHash(m_middleKmer, m_numHashes, m_k);
	RollingHash rightKmerHash(m_rightKmer, m_numHashes, m_k);

	ASSERT_EQ(rightKmerHash.getHash(), middleKmerHash.peekRight('A', 'C'));
}

TEST_F(RollingHashTest, peekLeft)
{
	RollingHash leftKmerHash(m_leftKmer, m_numHashes, m_k);
	RollingHash middleKmerHash(m_middleKmer, m_numHashes, m_k);

	ASSERT_EQ(leftKmerHash.getHash(), middleKmerHash.peekLeft('G', 'T'));
}

TEST_F(RollingHashTest, rollRight)
{
	RollingHash middleKmerHash(m_middleKmer, m_numHashes, m_k);
	RollingHash rightKmerHash(m_rightKmer, m_numHashes, m_k);

	middleKmerHash.rollRight('A', 'C');
	ASSERT_EQ(rightKmerHash.getHash(), middleKmerHash.getHash());
}

TEST_F(RollingHashTest, rollLeft)
{
	RollingHash leftKmerHash(m_leftKmer, m_numHashes, m_k);
	RollingHash middleKmerHash(m_middleKmer, m_numHashes, m_k);

	middleKmerHash.rollLeft('G', 'T');
	ASSERT_EQ(leftKmerHash.getHash(), middleKmerHash.getHash());
}

TEST_F(RollingHashTest, reset)
{
	RollingHash middleKmerHash(m_middleKmer, m_numHashes, m_k);
	RollingHash rightKmerHash(m_rightKmer, m_numHashes, m_k);

	middleKmerHash.reset(m_rightKmer);
	ASSERT_EQ(rightKmerHash.getHash(), middleKmerHash.getHash());
}
