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
