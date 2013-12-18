#include "connectpairs/BloomFilter.h"
#include "connectpairs/CountingBloomFilter.h"

#include <gtest/gtest.h>
#include <string>

using namespace std;

TEST(BloomFilter, base)
{
	BloomFilter x(100);
	EXPECT_EQ(x.size(), 100U);

	Kmer::setLength(16);
	Kmer a("AGATGTGCTGCCGCCT");
	Kmer b("TGGACAGCGTTACCTC");
	Kmer c("TAATAACAGTCCCTAT");
	Kmer d("GATCGTGGCGGGCGAT");

	x.insert(a);
	EXPECT_EQ(x.popcount(), 1U);
	EXPECT_TRUE(x[a]);
	x.insert(b);
	EXPECT_EQ(x.popcount(), 2U);
	EXPECT_TRUE(x[b]);
	EXPECT_TRUE(x[x.hash(b) % x.size()]);
	x.insert(x.hash(c) % x.size());
	EXPECT_EQ(x.popcount(), 3U);
	EXPECT_TRUE(x[c]);
	EXPECT_TRUE(x[x.hash(c) % x.size()]);

	EXPECT_FALSE(x[d]);
}

TEST(BloomFilter, serialization)
{
	BloomFilter origBloom(20);
	EXPECT_EQ(origBloom.size(), 20U);

	Kmer::setLength(16);
	Kmer a("AGATGTGCTGCCGCCT");
	Kmer b("TGGACAGCGTTACCTC");
	Kmer c("TAATAACAGTCCCTAT");

	origBloom.insert(a);
	origBloom.insert(b);
	origBloom.insert(c);

	EXPECT_TRUE(origBloom[a]);
	EXPECT_TRUE(origBloom[b]);
	EXPECT_TRUE(origBloom[c]);

	size_t origSize = origBloom.size();
	size_t origPopcount = origBloom.popcount();

	stringstream ss;
	ss << origBloom;
	ASSERT_TRUE(ss.good());

	BloomFilter copyBloom;
	ss >> copyBloom;
	ASSERT_TRUE(ss.good());

	EXPECT_EQ(copyBloom.size(), origSize);
	EXPECT_EQ(copyBloom.popcount(), origPopcount);

	EXPECT_TRUE(copyBloom[a]);
	EXPECT_TRUE(copyBloom[b]);
	EXPECT_TRUE(copyBloom[c]);
}

TEST(BloomFilter, loadAsUnion)
{
	size_t bits = 100;
	BloomFilter bloom1(bits);
	BloomFilter bloom2(bits);

	Kmer a("AGATGTGCTGCCGCCT");
	Kmer b("TGGACAGCGTTACCTC");

	bloom1.insert(a);
	bloom2.insert(b);

	EXPECT_TRUE(bloom1[a]);
	EXPECT_FALSE(bloom1[b]);
	EXPECT_FALSE(bloom2[a]);
	EXPECT_TRUE(bloom2[b]);

	BloomFilter unionBloom;
	
	stringstream ss;
	ss << bloom1;
	ASSERT_TRUE(ss.good());
	ss << bloom2;
	ASSERT_TRUE(ss.good());

	ss >> unionBloom;
	ss >> unionBloom;
	ASSERT_TRUE(ss.good());

	EXPECT_EQ(unionBloom.size(), bits);
	EXPECT_TRUE(unionBloom[a]);
	EXPECT_TRUE(unionBloom[b]);
}

TEST(CountingBloomFilter, base)
{
	CountingBloomFilter x(100);
	EXPECT_EQ(x.size(), 100U);

	Kmer::setLength(16);
	Kmer a("AGATGTGCTGCCGCCT");
	Kmer b("TGGACAGCGTTACCTC");
	Kmer c("TAATAACAGTCCCTAT");
	Kmer d("GATCGTGGCGGGCGAT");

	x.insert(a);
	EXPECT_EQ(x.popcount(), 0U);
	EXPECT_FALSE(x[a]);
	x.insert(a);
	EXPECT_EQ(x.popcount(), 1U);
	EXPECT_TRUE(x[a]);
	x.insert(b);
	EXPECT_EQ(x.popcount(), 1U);
	EXPECT_FALSE(x[b]);
	x.insert(c);
	EXPECT_EQ(x.popcount(), 1U);
	EXPECT_FALSE(x[c]);
	x.insert(b);
	EXPECT_EQ(x.popcount(), 2U);
	EXPECT_TRUE(x[b]);
	EXPECT_TRUE(x[x.hash(b) % x.size()]);
	x.insert(x.hash(c) % x.size());
	EXPECT_EQ(x.popcount(), 3U);
	EXPECT_TRUE(x[c]);

	EXPECT_FALSE(x[d]);
}
