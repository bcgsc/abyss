#include "PairedDBG/KmerPair.h"

#include <gtest/gtest.h>

using namespace std;

string seq1 = "AACCTTGG";
string seq2 = "ACGTACGT";
string seq = "AACCTTGGNNNNNACGTACGT";

TEST(KmerPair, constructors)
{
	Kmer::setLength(8);
	KmerPair::setLength(21);

	Kmer kmer1(seq1);
	Kmer kmer2(seq2);
	KmerPair k1(kmer1, kmer2);
	KmerPair k2(seq1, seq2);
	KmerPair k3(seq);
	KmerPair k4(seq1, seq1);

	EXPECT_EQ(k1, k2);
	EXPECT_EQ(k1, k3);
	EXPECT_NE(k1, k4);

	EXPECT_EQ(KmerPair::length(), 21u);
}

TEST(KmerPair, str)
{
	Kmer::setLength(8);
	KmerPair::setLength(21);

	KmerPair k(seq);
	string res1 = k.str();
	EXPECT_EQ(res1, seq);

	KmerPair::setLength(22); // Maybe this shouldn't be allowed?
	string res2 = k.str();
	EXPECT_EQ(res2, "AACCTTGGNNNNNNACGTACGT"); // Add one 'N' from seq
}

TEST(KmerPair, reverseComplement)
{
	string rcseq1 = "CCAAGGTT";
	string rcseq2 = "ACGTACGT";
	EXPECT_EQ(rcseq1, reverseComplement(seq1)); // just to make sure ;)

	KmerPair k(seq1, seq2);
	KmerPair rck(rcseq2, rcseq1);
	EXPECT_EQ(rck, reverseComplement(k));
	k.reverseComplement();
	EXPECT_EQ(rck, k);
}

TEST(KmerPair, isPalindrome)
{
	Kmer::setLength(8);
	KmerPair::setLength(21);

	string rcseq1 = reverseComplement(seq1);
	KmerPair kp(seq1, rcseq1);
	EXPECT_EQ(kp, reverseComplement(kp));
	EXPECT_TRUE(kp.isPalindrome());

	string pal("AGAATTCT");
	Kmer k(pal);
	EXPECT_TRUE(k.isPalindrome());
	KmerPair kp_pal(k, k);
	EXPECT_TRUE(k.isPalindrome());

	KmerPair kp_npal(pal, seq2);
	EXPECT_FALSE(kp_npal.isPalindrome());
}

TEST(KmerPair, isPalindrome_edge)
{
	Kmer::setLength(4);

	string epal = "CCGCNNNNAGCG";
	KmerPair kp(epal);
	EXPECT_FALSE(kp.isPalindrome());
	EXPECT_FALSE(kp.isPalindrome(ANTISENSE));
	EXPECT_TRUE(kp.isPalindrome(SENSE));
}

