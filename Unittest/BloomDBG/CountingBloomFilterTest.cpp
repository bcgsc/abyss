#include "vendor/btl_bloomfilter/CountingBloomFilter.hpp"
#include "BloomDBG/RollingHashIterator.h"

#include <gtest/gtest.h>

using namespace std;
typedef uint64_t hash_t;

TEST(CountingBloomFilter, base) {
  const unsigned bloomSize = 1000;
  const unsigned numHashes = 1;
  const unsigned numLevels = 2;
  const unsigned k = 16;

  CountingBloomFilter<uint8_t> x(bloomSize, numHashes, k, numLevels);
  EXPECT_EQ(x.sizeInBytes(), bloomSize);

  const char *a = "AGATGTGCTGCCGCCT";
  const char *b = "TGGACAGCGTTACCTC";
  const char *c = "TAATAACAGTCCCTAT";
  const char *d = "GATCGTGGCGGGCGAT";
  const char *e = "TTTTTTTTTTTTTTTT";

  RollingHashIterator itA(a, numHashes, k);
  RollingHashIterator itB(b, numHashes, k);
  RollingHashIterator itC(c, numHashes, k);
  RollingHashIterator itD(d, numHashes, k);
  RollingHashIterator itE(e, numHashes, k);

  x.insert(*itA);
  EXPECT_EQ(x.filtered_popcount(), 0U);
  EXPECT_FALSE(x.contains(*itE));
  x.insert(*itA);
  EXPECT_EQ(x.filtered_popcount(), 1U);
  EXPECT_TRUE(x.contains(*itA));
  x.insert(*itB);
  EXPECT_EQ(x.filtered_popcount(), 1U);
  EXPECT_FALSE(x.contains(*itB));
  x.insert(*itC);
  EXPECT_EQ(x.filtered_popcount(), 1U);
  EXPECT_FALSE(x.contains(*itC));
  x.insert(*itB);
  EXPECT_EQ(x.filtered_popcount(), 2U);
  EXPECT_TRUE(x.contains(*itB));
  EXPECT_FALSE(x.contains(*itD));
}
