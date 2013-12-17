#include "connectpairs/DBGBloomAlgorithms.h"
#include "Common/Sequence.h"
#include <gtest/gtest.h>
#include <string>

using namespace std;

namespace {

/**
 * Test fixture for the getStartKmerPos function.
 *
 * getStartKmerPos chooses the kmer at the
 * end of the longest string of kmer matches with
 * within the read.
 *
 * Once the longest match region is identified,
 * the kmer should be selected at the end of the
 * region that is closest to the gap between
 * read pairs.
 */
class GetStartKmerPosTest : public ::testing::Test {

protected:

	FastaRecord testRead;
	static const int k = 2;
	// set this large to avoid false positives
	static const int bloomFilterSize = 1000000;

	GetStartKmerPosTest() {
		Kmer::setLength(k);
		testRead.seq = "TACAGTG";
	}

};

TEST_F(GetStartKmerPosTest, FullReadMatch)
{
	BloomFilter bloom(bloomFilterSize);
	bloom.loadSeq(k, testRead.seq);
	DBGBloom g(bloom);

	EXPECT_EQ(testRead.seq.length() - k, getStartKmerPos(k, testRead, g));
	// true indicates revese complement (second read)
	EXPECT_EQ(getStartKmerPos(k, testRead, g), getStartKmerPos(k, testRead, g, true));
}

TEST_F(GetStartKmerPosTest, FullReadMismatch)
{
	BloomFilter bloom(bloomFilterSize);
	DBGBloom g(bloom);
	EXPECT_EQ(NO_MATCH, getStartKmerPos(k, testRead, g));
}

TEST_F(GetStartKmerPosTest, SelectLongestMatchRegion)
{
	const string& seq = testRead.seq;
	BloomFilter bloom(bloomFilterSize);
	DBGBloom g(bloom);

	// This loop creates kmer match vector 101101
	for (unsigned i = 0; i < seq.length() - k + 1; i++) {
		// non-matching kmers
		if (i == 1 || i == 4)
			continue;
		bloom.loadSeq(k, seq.substr(i,k));
	}

	EXPECT_EQ(3u, getStartKmerPos(k, testRead, g));
	// true indicates revese complement (second read)
	EXPECT_EQ(getStartKmerPos(k, testRead, g),
		getStartKmerPos(k, testRead, g, true));
}

TEST_F(GetStartKmerPosTest, EqualLengthMatchRegions)
{
	const string& seq = testRead.seq;
	BloomFilter bloom(bloomFilterSize);
	DBGBloom g(bloom);

	// This loop creates kmer match vector 110011
	for (unsigned i = 0; i < seq.length() - k + 1; i++) {
		// non-matching kmers
		if (i == 0 || i == 3)
			continue;
		bloom.loadSeq(k, seq.substr(i,k));
	}

	EXPECT_EQ(5u, getStartKmerPos(k, testRead, g));
	// true indicates revese complement (second read)
	EXPECT_EQ(getStartKmerPos(k, testRead, g),
		getStartKmerPos(k, testRead, g, true));
}

}

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

	BloomFilter bloom(1000);
	DBGBloom g(bloom);

	bloom.loadSeq(k, read1.seq);
	bloom.loadSeq(k, read2.seq);

	vector<FastaRecord> mergedSeqs;
	SearchResult result = connectPairs(k, read1, read2, g, 1, 0, 4);

	EXPECT_EQ(FOUND_PATH, result.pathResult);
	ASSERT_EQ(1u, result.mergedSeqs.size());
	EXPECT_EQ("GATG", result.mergedSeqs[0].seq);
}
