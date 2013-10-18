#include "Common/Sequence.h"

#include <gtest/gtest.h>

using namespace std;

TEST(reverseComplement, base)
{
	string s = "AGATGTGCTGCCGCCTTGGACAGCGTTACCTCTAATAACAGTCCCTATGA";
	string rc = "TCATAGGGACTGTTATTAGAGGTAACGCTGTCCAAGGCGGCAGCACATCT";
	EXPECT_EQ(reverseComplement(s), rc);
	EXPECT_EQ(reverseComplement(reverseComplement(s)), s);
}
