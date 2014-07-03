#include "Konnector/konnector.h"
#include <iostream>

#include <gtest/gtest.h>

using namespace std;

// workaround: opt::k must be defined because
// it is used by write_dot(..)
namespace opt { unsigned k = 0; }

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

TEST(ConnectPairsTest, MergeOverlappingPair)
{
	// Merged seq: GATG
	// Read 1:     GAT  =>
	// Read 2:      TAC <=

	const int k = 2;
	Kmer::setLength(k);
	const int readLength = 3;
	string mergedSeq = "GATG";

	FastaRecord read1, read2;
	read1.id = "read/1";
	read1.seq = mergedSeq.substr(0,readLength);
	read2.id = "read/2";
	read2.seq = reverseComplement(mergedSeq.substr(1,readLength));

	BloomFilter bloom(1000);
	DBGBloom g(bloom);

	bloom.loadSeq(k, read1.seq);
	bloom.loadSeq(k, read2.seq);

	vector<FastaRecord> mergedSeqs;

	ConnectPairsParams params;
	params.maxPaths = 1;
	params.minMergedSeqLen = 0;
	params.maxMergedSeqLen = 4;

	ConnectPairsResult result = connectPairs(k, read1, read2, g, params);

	EXPECT_EQ(FOUND_PATH, result.pathResult);
	ASSERT_EQ(1u, result.mergedSeqs.size());
	EXPECT_EQ("GATG", result.mergedSeqs[0].seq);

}
