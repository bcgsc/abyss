#include "Bloom/Bloom.h"
#include "Bloom/KonnectorBloomFilter.h"
#include "Bloom/CascadingBloomFilter.h"
#include "BloomDBG/RollingHashIterator.h"
#include "Common/BitUtil.h"

#include <gtest/gtest.h>
#include <string>

using namespace std;

TEST(BloomFilter, base)
{
	KonnectorBloomFilter x(10000, 16);
	EXPECT_EQ(x.size(), 10000U);

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
	RollingHashIterator itb("TGGACAGCGTTACCTC", 1, 16);
	EXPECT_TRUE(x[*itb]);
	RollingHashIterator itc("TAATAACAGTCCCTAT", 1, 16);
	x.insert(*itc);
	EXPECT_EQ(x.popcount(), 3U);
	EXPECT_TRUE(x[c]);

	EXPECT_FALSE(x[d]);
}

TEST(BloomFilter, serialization)
{
	KonnectorBloomFilter origBloom(24, 16);
	EXPECT_EQ(origBloom.size(), 24U);

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

	KonnectorBloomFilter copyBloom;
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
	size_t bits = 10000;
	KonnectorBloomFilter bloom1(bits, 16);
	KonnectorBloomFilter bloom2(bits, 16);

	Kmer a("AGATGTGCTGCCGCCT");
	Kmer b("TGGACAGCGTTACCTC");

	bloom1.insert(a);
	bloom2.insert(b);

	EXPECT_TRUE(bloom1[a]);
	EXPECT_FALSE(bloom1[b]);
	EXPECT_FALSE(bloom2[a]);
	EXPECT_TRUE(bloom2[b]);

	KonnectorBloomFilter unionBloom;

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
	size_t bits = 10000;
	KonnectorBloomFilter bloom1(bits, 16);
	KonnectorBloomFilter bloom2(bits, 16);

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

	KonnectorBloomFilter intersectBloom;

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
	CascadingBloomFilter x(10000, 2, 16);
	EXPECT_EQ(x.size(), 10000U);

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
	RollingHashIterator itb("TGGACAGCGTTACCTC", 1, 16);
	EXPECT_TRUE(x[(*itb)[0]]);
	RollingHashIterator itc("TAATAACAGTCCCTAT", 1, 16);
	x.insert(*itc);
	EXPECT_EQ(x.popcount(), 3U);
	EXPECT_TRUE(x[c]);

	EXPECT_FALSE(x[d]);
}

