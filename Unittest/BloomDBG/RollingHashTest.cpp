#include "BloomDBG/RollingHash.h"

#include <gtest/gtest.h>

using namespace std;

/** test fixture for RollingHash tests */
class RollingHashTest : public ::testing::Test
{
protected:

	const unsigned m_numHashes;
	const unsigned m_k;
	const string m_kmer;
	const string m_nextKmer;

	RollingHashTest() : m_numHashes(2), m_k(4), m_kmer("ACGT"),
		m_nextKmer("CGTC") {}
};

TEST_F(RollingHashTest, peekRight)
{
	RollingHash kmerHash(m_kmer, m_numHashes, m_k);
	RollingHash nextKmerHash(m_nextKmer, m_numHashes, m_k);

	ASSERT_EQ(nextKmerHash.getHash(), kmerHash.peekRight('A', 'C'));
}

TEST_F(RollingHashTest, rollRight)
{
	RollingHash kmerHash(m_kmer, m_numHashes, m_k);
	RollingHash nextKmerHash(m_nextKmer, m_numHashes, m_k);

	kmerHash.rollRight('A', 'C');
	ASSERT_EQ(nextKmerHash.getHash(), kmerHash.getHash());
}

TEST_F(RollingHashTest, reset)
{
	RollingHash kmerHash(m_kmer, m_numHashes, m_k);
	RollingHash nextKmerHash(m_nextKmer, m_numHashes, m_k);

	kmerHash.reset(m_nextKmer);
	ASSERT_EQ(nextKmerHash.getHash(), kmerHash.getHash());
}
