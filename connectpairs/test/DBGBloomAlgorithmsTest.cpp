#include "connectpairs/DBGBloomAlgorithms.h"
#include "Common/Sequence.h"
#include <gtest/gtest.h>
#include <string>

using namespace std;

TEST(DBGBloomAlgorithmsTest, MergeOverlappingPair)
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

	DBGBloom g(k, mergedSeq.size());
	g.assign(read1.seq);
	g.assign(read2.seq);

	vector<FastaRecord> mergedSeqs;
	PathSearchResult result;
	result = connectPairs(read1, read2, g, mergedSeqs, 1, 4);

	EXPECT_EQ(FOUND_PATH, result);
	ASSERT_EQ(1, mergedSeqs.size());
	EXPECT_EQ("GATG", mergedSeqs[0].seq);
}
