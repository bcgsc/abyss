#include "BloomDBG/MaskedKmer.h"

#include <gtest/gtest.h>

using namespace std;

TEST(MaskedKmerTest, trivialMask)
{
	Kmer::setLength(4);

	MaskedKmer kmer1("ACGT");
	MaskedKmer kmer2("ACGT");

	ASSERT_EQ(kmer1, kmer2);
}

TEST(MaskedKmerTest, nonTrivialMask)
{
	Kmer::setLength(4);

	const string mask("1001");
	MaskedKmer kmer1("ACGT", mask);
	MaskedKmer kmer2("ATTT", mask);

	ASSERT_EQ(kmer1, kmer2);
}
