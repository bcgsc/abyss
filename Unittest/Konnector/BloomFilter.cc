#include "Bloom/Bloom.h"
#include "Bloom/BloomFilter.h"
#include "Bloom/CascadingBloomFilter.h"
#include "Bloom/BloomFilterWindow.h"
#include "Bloom/CascadingBloomFilterWindow.h"
#include "Common/BitUtil.h"

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
	EXPECT_TRUE(x[Bloom::hash(b) % x.size()]);
	x.insert(Bloom::hash(c) % x.size());
	EXPECT_EQ(x.popcount(), 3U);
	EXPECT_TRUE(x[c]);
	EXPECT_TRUE(x[Bloom::hash(c) % x.size()]);

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
	unionBloom.read(ss, BITWISE_OR);
	ASSERT_TRUE(ss.good());

	EXPECT_EQ(unionBloom.size(), bits);
	EXPECT_TRUE(unionBloom[a]);
	EXPECT_TRUE(unionBloom[b]);
}

TEST(BloomFilter, intersect)
{
	size_t bits = 100;
	BloomFilter bloom1(bits);
	BloomFilter bloom2(bits);

	Kmer a("AGATGTGCTGCCGCCT");
	Kmer b("TGGACAGCGTTACCTC");
	Kmer c("AGCTAGCTAGCTAGCT");

	bloom1.insert(a);
	bloom2.insert(b);

	bloom1.insert(c);
	bloom2.insert(c);

	EXPECT_TRUE(bloom1[a]);
	EXPECT_TRUE(bloom1[c]);
	EXPECT_FALSE(bloom1[b]);
	EXPECT_FALSE(bloom2[a]);
	EXPECT_TRUE(bloom2[b]);
	EXPECT_TRUE(bloom2[c]);

	BloomFilter intersectBloom;

	stringstream ss;
	ss << bloom1;
	ASSERT_TRUE(ss.good());
	ss << bloom2;
	ASSERT_TRUE(ss.good());

	ss >> intersectBloom;
	ASSERT_TRUE(ss.good());
	intersectBloom.read(ss, BITWISE_AND);
	ASSERT_TRUE(ss.good());

	EXPECT_EQ(intersectBloom.size(), bits);
	EXPECT_FALSE(intersectBloom[a]);
	EXPECT_FALSE(intersectBloom[b]);
	EXPECT_TRUE(intersectBloom[c]);
}

TEST(CascadingBloomFilter, base)
{
	CascadingBloomFilter x(100);
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
	EXPECT_TRUE(x[Bloom::hash(b) % x.size()]);
	x.insert(Bloom::hash(c) % x.size());
	EXPECT_EQ(x.popcount(), 3U);
	EXPECT_TRUE(x[c]);

	EXPECT_FALSE(x[d]);
}

TEST(BloomFilter, windowSerialization)
{
	size_t bits = 100;
	size_t pos = 80;

	BloomFilterWindow window(bits, bits/2, bits-1);

	window.insert(pos);
	EXPECT_TRUE(window[pos]);

	stringstream ss;
	ss << window;
	ASSERT_TRUE(ss.good());

	BloomFilterWindow windowCopy;
	ss >> windowCopy;
	ASSERT_TRUE(ss.good());
	EXPECT_TRUE(windowCopy[pos]);
	EXPECT_EQ(window.fullBloomSize(), windowCopy.fullBloomSize());
	EXPECT_EQ(window.startBitPos(), windowCopy.startBitPos());
	EXPECT_EQ(window.endBitPos(), windowCopy.endBitPos());
}

TEST(BloomFilter, windowUnion)
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
	unionBloom.read(ss, BITWISE_OR);
	ASSERT_TRUE(ss.good());

	EXPECT_EQ(2U, unionBloom.popcount());
	EXPECT_TRUE(unionBloom[pos1]);
	EXPECT_TRUE(unionBloom[pos2]);
}

TEST(CascadingBloomFilter, window)
{
	size_t bits = 100;
	size_t pos1 = 25;
	size_t pos2 = 80;
	size_t pos3 = 50;

	CascadingBloomFilter CascadingBloom(bits);

	// set a bit in both halves of the second level
	// bloom filter
	CascadingBloom.insert(pos1);
	CascadingBloom.insert(pos1);
	CascadingBloom.insert(pos2);
	CascadingBloom.insert(pos2);
	CascadingBloom.insert(pos3);

	EXPECT_TRUE(CascadingBloom[pos1]);
	EXPECT_TRUE(CascadingBloom[pos2]);
	EXPECT_EQ(2U, CascadingBloom.getBloomFilter(1).popcount());

	CascadingBloomFilterWindow window1(bits, 0, bits/2 - 1);
	window1.insert(pos1);
	window1.insert(pos1);

	CascadingBloomFilterWindow window2(bits, bits/2, bits - 1);
	window2.insert(pos2);
	window2.insert(pos2);

	window2.insert(pos3);
	window1.insert(pos3);

	stringstream ss;
	ss << window1;
	ASSERT_TRUE(ss.good());
	ss << window2;
	ASSERT_TRUE(ss.good());

	BloomFilter unionBloom;
	ss >> unionBloom;
	ASSERT_TRUE(ss.good());
	unionBloom.read(ss, BITWISE_OR);
	ASSERT_TRUE(ss.good());

	EXPECT_EQ(2U, unionBloom.popcount());
	EXPECT_TRUE(unionBloom[pos1]);
	EXPECT_TRUE(unionBloom[pos2]);
	EXPECT_FALSE(unionBloom[pos3]);
}
