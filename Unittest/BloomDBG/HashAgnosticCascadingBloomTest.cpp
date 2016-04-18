#include "BloomDBG/RollingHashIterator.h"
#include "BloomDBG/HashAgnosticCascadingBloom.h"

#include <gtest/gtest.h>

using namespace std;
typedef uint64_t hash_t;

TEST(HashAgnosticCascadingBloom, base)
{
	const unsigned bloomSize = 1000;
	const unsigned numHashes = 1;
	const unsigned numLevels = 2;
	const unsigned k = 16;

	HashAgnosticCascadingBloom x(bloomSize, numHashes, numLevels, k);
	EXPECT_EQ(x.size(), bloomSize);

	const char* a = "AGATGTGCTGCCGCCT";
	const char* b = "TGGACAGCGTTACCTC";
	const char* c = "TAATAACAGTCCCTAT";
	const char* d = "GATCGTGGCGGGCGAT";

	RollingHashIterator itA(a, numHashes, k);
	RollingHashIterator itB(b, numHashes, k);
	RollingHashIterator itC(c, numHashes, k);
	RollingHashIterator itD(d, numHashes, k);
	size_t hash;

	x.insert(*itA);
	EXPECT_EQ(x.popcount(), 0U);
	EXPECT_FALSE(x.contains(&hash));
	x.insert(*itA);
	EXPECT_EQ(x.popcount(), 1U);
	EXPECT_TRUE(x.contains(*itA));
	x.insert(*itB);
	EXPECT_EQ(x.popcount(), 1U);
	EXPECT_FALSE(x.contains(*itB));
	x.insert(*itC);
	EXPECT_EQ(x.popcount(), 1U);
	EXPECT_FALSE(x.contains(*itC));
	x.insert(*itB);
	EXPECT_EQ(x.popcount(), 2U);
	EXPECT_TRUE(x.contains(*itB));
	EXPECT_FALSE(x.contains(*itD));
}
