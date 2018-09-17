#include "BloomDBG/LightweightKmer.h"
#include "Common/Kmer.h"

#include <gtest/gtest.h>

TEST(LightweightKmerTest, isCanonical)
{
	Kmer::setLength(5);

	const LightweightKmer kmer1("ACGTA");
	const LightweightKmer kmer2("TACGT");

	EXPECT_TRUE(kmer1.isCanonical());
	EXPECT_FALSE(kmer2.isCanonical());
}

TEST(LightweightKmerTest, canonicalize)
{
	Kmer::setLength(5);

	const LightweightKmer kmer1("ACGTA");
	const LightweightKmer kmer2("TACGT");

	LightweightKmer kmer1Copy(kmer1.c_str());
	kmer1Copy.canonicalize();
	EXPECT_EQ(kmer1, kmer1Copy);

	LightweightKmer kmer2Copy(kmer2.c_str());
	kmer2Copy.canonicalize();
	EXPECT_NE(kmer2, kmer2Copy);
	EXPECT_EQ(kmer1, kmer2Copy);
}

TEST(LightweightKmerTest, LessThan)
{
	/*
	 * Note: the less-than operator (i.e. `operator<`)
	 * is written to be invariant under reverse-complement.
	 * In other words, it compares the canonical orientations
	 * of the two k-mers.
	 */

	Kmer::setLength(5);

	const LightweightKmer kmer1("ACGTA");
	const LightweightKmer rcKmer1("TACGT");

	const LightweightKmer kmer2("TGCAT");
	const LightweightKmer rcKmer2("ATGCA");

	ASSERT_TRUE(kmer1.isCanonical());
	ASSERT_FALSE(kmer2.isCanonical());

	/*
	 * we expect kmer1 < kmer2, regardless of the orientation
	 * of the two k-mers
	 */

	ASSERT_TRUE(kmer1 < kmer2);
	ASSERT_TRUE(kmer1 < rcKmer2);
	ASSERT_TRUE(rcKmer1 < kmer2);
	ASSERT_TRUE(rcKmer1 < rcKmer2);
}
