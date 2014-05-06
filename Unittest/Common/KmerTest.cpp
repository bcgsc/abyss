#include "Common/Kmer.h"
#include <gtest/gtest.h>
#include <iostream>

TEST(Kmer, canonicalize)
{
	Kmer::setLength(4);
	Kmer canonical("ATGC");
	Kmer nonCanonical("GCAT");
	Kmer palindrome("ACGT");

	Kmer kmer = canonical;
	kmer.canonicalize();
	EXPECT_EQ(canonical, kmer);

	kmer = nonCanonical;
	kmer.canonicalize();
	EXPECT_EQ(canonical, kmer);

	kmer = palindrome;
	kmer.canonicalize();
	EXPECT_EQ(palindrome, kmer);

	Kmer::setLength(5);
	Kmer oddLength("GCTCG");
	Kmer oddLengthCanonical("CGAGC");

	kmer = oddLength;
	kmer.canonicalize();
	EXPECT_EQ(oddLengthCanonical, kmer);
}

