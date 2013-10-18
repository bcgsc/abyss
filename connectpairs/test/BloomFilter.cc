#include "connectpairs/BloomFilter.h"

#include <gtest/gtest.h>
#include <string>

using namespace std;

TEST(BloomFilter, base)
{
	BloomFilter x(100);
	EXPECT_EQ(x.size(), 100);

	Kmer::setLength(16);
	Kmer a("AGATGTGCTGCCGCCT");
	Kmer b("TGGACAGCGTTACCTC");
	Kmer c("TAATAACAGTCCCTAT");
	Kmer d("GATCGTGGCGGGCGAT");

	x.insert(a);
	EXPECT_EQ(x.popcount(), 1);
	EXPECT_TRUE(x[a]);
	x.insert(b);
	EXPECT_EQ(x.popcount(), 2);
	EXPECT_TRUE(x[b]);
	x.insert(c);
	EXPECT_EQ(x.popcount(), 3);
	EXPECT_TRUE(x[c]);
	EXPECT_TRUE(x[x.hash(c) % x.size()]);

	EXPECT_FALSE(x[d]);
}
