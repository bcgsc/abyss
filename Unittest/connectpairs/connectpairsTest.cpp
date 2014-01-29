#include "connectpairs/connectpairs.h"

#include <gtest/gtest.h>

TEST(maskNew, read1)
{
	FastqRecord r1("1", "", "ACGTACGT", "BBBBBBBB");
	FastqRecord r2;
	FastaRecord read("2", "", "ACGTACGT");

	EXPECT_TRUE(maskNew(r1, r2, read) == 0u);
	EXPECT_TRUE(read.seq == "ACGTACGT");

	read = FastaRecord("2", "", "ACGTACGTA");
	EXPECT_TRUE(maskNew(r1, r2, read);
	EXPECT_TRUE(read.seq == "ACGTACGTa");
}
