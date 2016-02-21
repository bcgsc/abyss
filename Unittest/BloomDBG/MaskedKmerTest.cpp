#include "BloomDBG/MaskedKmer.h"

#include <gtest/gtest.h>

using namespace std;

TEST(MaskedKmerTest, trivialMask)
{
	MaskedKmer::setLength(4);

	MaskedKmer kmer1("ACGT");
	MaskedKmer kmer2("ACGT");

	ASSERT_EQ(kmer1, kmer2);
}

TEST(MaskedKmerTest, nonTrivialMask)
{
	MaskedKmer::setLength(4);
	MaskedKmer::setMask("1001");

	MaskedKmer kmer1("ACGT");
	MaskedKmer kmer2("ATTT");

	ASSERT_EQ(kmer1, kmer2);
}
