#include "connectpairs/BloomFilter.h"
#include "connectpairs/CountingBloomFilter.h"
#include "connectpairs/BloomFilterWindow.h"
#include "connectpairs/CountingBloomFilterWindow.h"

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

	EXPECT_EQ(origSize, copyBloom.size());
	EXPECT_EQ(origPopcount, copyBloom.popcount());

	EXPECT_TRUE(copyBloom[a]);
	EXPECT_TRUE(copyBloom[b]);
	EXPECT_TRUE(copyBloom[c]);
}

TEST(BloomFilter, union_)
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
	ASSERT_TRUE(ss.good());
	unionBloom.read(ss, true);
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

TEST(BloomFilter, shrink)
{
	BloomFilter big(10);
	BloomFilter small;

	big.insert(1);
	big.insert(8);

	EXPECT_EQ(2U, big.popcount());
	EXPECT_TRUE(big[1]);
	EXPECT_TRUE(big[8]);

	stringstream ss;
	ss << big;
	ASSERT_TRUE(ss.good());
	small.read(ss, false, 2);
	ASSERT_TRUE(ss.good());

	EXPECT_EQ(5U, small.size());
	EXPECT_EQ(2U, small.popcount());
	EXPECT_TRUE(small[1]);
	EXPECT_TRUE(small[3]);
}

TEST(BloomFilter, window)
{
	size_t bits = 100;
	size_t pos1 = 25;
	size_t pos2 = 80;

	BloomFilter bloom(bits);

	// set a bit in both halves of bloom filter
	// bit array
	bloom.insert(pos1);
	bloom.insert(pos2);

	EXPECT_TRUE(bloom[pos1]);
	EXPECT_TRUE(bloom[pos2]);
	EXPECT_EQ(2U, bloom.popcount());

	BloomFilterWindow window1(bits, 0, bits/2 - 1);
	window1.insert(pos1);
	EXPECT_TRUE(window1[pos1]);

	BloomFilterWindow window2(bits, bits/2, bits - 1);
	window2.insert(pos2);
	EXPECT_TRUE(window2[pos2]);

	stringstream ss;
	ss << window1;
	ASSERT_TRUE(ss.good());
	ss << window2;
	ASSERT_TRUE(ss.good());

	BloomFilter unionBloom;
	ss >> unionBloom;
	ASSERT_TRUE(ss.good());
	// true means load the union
	unionBloom.read(ss, true);
	ASSERT_TRUE(ss.good());

	EXPECT_EQ(2U, unionBloom.popcount());
	EXPECT_TRUE(unionBloom[pos1]);
	EXPECT_TRUE(unionBloom[pos2]);
}

TEST(CountingBloomFilter, window)
{
	size_t bits = 100;
	size_t pos1 = 25;
	size_t pos2 = 80;

	CountingBloomFilter countingBloom(bits);

	// set a bit in both halves of the second level
	// bloom filter
	countingBloom.insert(pos1);
	countingBloom.insert(pos1);
	countingBloom.insert(pos2);
	countingBloom.insert(pos2);

	EXPECT_TRUE(countingBloom[pos1]);
	EXPECT_TRUE(countingBloom[pos2]);
	EXPECT_EQ(2U, countingBloom.getBloomFilter(1).popcount());

	CountingBloomFilterWindow window1(bits, 0, bits/2 - 1);
	window1.insert(pos1);
	window1.insert(pos1);

	CountingBloomFilterWindow window2(bits, bits/2, bits - 1);
	window2.insert(pos2);
	window2.insert(pos2);

	stringstream ss;
	ss << window1.getBloomFilter(1);
	ASSERT_TRUE(ss.good());
	ss << window2.getBloomFilter(1);
	ASSERT_TRUE(ss.good());

	BloomFilter unionBloom;
	ss >> unionBloom;
	ASSERT_TRUE(ss.good());
	// true means load the union
	unionBloom.read(ss, true);
	ASSERT_TRUE(ss.good());

	EXPECT_EQ(2U, unionBloom.popcount());
	EXPECT_TRUE(unionBloom[pos1]);
	EXPECT_TRUE(unionBloom[pos2]);
}
