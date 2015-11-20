#include "BloomDBG/RollingHashIterator.h"
#include "BloomDBG/HashAgnosticCascadingBloom.h"
#include "lib/bloomfilter-521e80c5c619a9a8e3d6389dc3b597a75bdf2aaa/BloomFilter.hpp"

#include <gtest/gtest.h>

using namespace std;
typedef uint64_t hash_t;

TEST(HashAgnosticCascadingBloom, base)
{
	const unsigned bloomSize = 1000;
	const unsigned numHashes = 2;
	const unsigned numLevels = 2;
	const unsigned k = 16;

	HashAgnosticCascadingBloom x(bloomSize, numHashes, numLevels, k);
	EXPECT_EQ(x.size(), bloomSize);

	const char* a = "AGATGTGCTGCCGCCT";
	const char* b = "TGGACAGCGTTACCTC";
	const char* c = "TAATAACAGTCCCTAT";
	const char* d = "GATCGTGGCGGGCGAT";

	RollingHashIterator itA(a, k, numHashes);
	RollingHashIterator itB(b, k, numHashes);
	RollingHashIterator itC(c, k, numHashes);
	RollingHashIterator itD(d, k, numHashes);

	x.insert(*itA);
	EXPECT_EQ(x.popcount(), 0U);
	EXPECT_FALSE(x.contains(*itA));
	x.insert(*itA);
	EXPECT_EQ(x.popcount(), 2U);
	EXPECT_TRUE(x.contains(*itA));
	x.insert(*itB);
	EXPECT_EQ(x.popcount(), 2U);
	EXPECT_FALSE(x.contains(*itB));
	x.insert(*itC);
	EXPECT_EQ(x.popcount(), 2U);
	EXPECT_FALSE(x.contains(*itC));
	x.insert(*itB);
	EXPECT_EQ(x.popcount(), 4U);
	EXPECT_TRUE(x.contains(*itB));
	EXPECT_FALSE(x.contains(*itD));
}
