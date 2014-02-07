#include "connectpairs/connectpairs.h"
#include <iostream>

#include <gtest/gtest.h>

using namespace std;

TEST(maskNew, read1)
{
	FastqRecord r1("1", "", "ACGTACGT", "BBBBBBBB");
	FastqRecord r2;
	FastaRecord read("2", "", "ACGTACGT");

	int mask = 1;
	EXPECT_TRUE(maskNew(r1, r2, read, mask) == 0u);
	EXPECT_TRUE(read.seq == "ACGTACGT");

	read = FastaRecord("2", "", "ACGTACGTA");
	EXPECT_TRUE(maskNew(r1, r2, read, mask) == 0u);
	cout << read.seq << endl;
	EXPECT_TRUE(read.seq == "ACGTACGTa");
}

TEST(maskNew, mask)
{
	FastqRecord r1("1", "", "ACGTA", "BBBBB");
	FastqRecord r2;
	FastaRecord read("2", "", "ACGTACGT");

	int mask = 0;
	EXPECT_TRUE(maskNew(r1, r2, read, mask) == 0u);
	EXPECT_TRUE(read.seq == "ACGTACGT");
}
