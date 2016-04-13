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

	RollingHashIterator itA(a, k);
	RollingHashIterator itB(b, k);
	RollingHashIterator itC(c, k);
	RollingHashIterator itD(d, k);
	size_t hash;

	hash = *itA;
	x.insert(&hash);
	EXPECT_EQ(x.popcount(), 0U);
	EXPECT_FALSE(x.contains(&hash));
	x.insert(&hash);
	EXPECT_EQ(x.popcount(), 1U);
	EXPECT_TRUE(x.contains(&hash));
	hash = *itB;
	x.insert(&hash);
	EXPECT_EQ(x.popcount(), 1U);
	EXPECT_FALSE(x.contains(&hash));
	hash = *itC;
	x.insert(&hash);
	EXPECT_EQ(x.popcount(), 1U);
	EXPECT_FALSE(x.contains(&hash));
	hash = *itB;
	x.insert(&hash);
	EXPECT_EQ(x.popcount(), 2U);
	EXPECT_TRUE(x.contains(&hash));
	hash = *itD;
	EXPECT_FALSE(x.contains(&hash));
}
