#include "connectpairs/DBGBloomAlgorithms.h"
#include "Common/Sequence.h"
#include <gtest/gtest.h>
#include <string>

using namespace std;

/**
 * Test fixture for the getStartKmerPos function.
 *
 * getStartKmerPos chooses the kmer at the
 * beginning of the longest series kmer matches with
 * within the read.  In the case where there are
 * two equal-length series of matches,
 * getStartKmerPos choose the kmer position closest
 * to the beginning (5' end) of the read.
 */
class GetStartKmerPosTest : public testing::Test {

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

	EXPECT_EQ(0U, getStartKmerPos(k, testRead, g, false, true));
	// true indicates revese complement (second read)
	EXPECT_EQ(0U, getStartKmerPos(k, testRead, g, true, true));
}

TEST_F(GetStartKmerPosTest, FullReadMismatch)
{
	BloomFilter bloom(bloomFilterSize);
	// Leave the bloom filter empty to generate a mismatch
	// for every kmer in the read.
	DBGBloom g(bloom);
	EXPECT_EQ(NO_MATCH, getStartKmerPos(k, testRead, g, false, true));
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

	EXPECT_EQ(2U, getStartKmerPos(k, testRead, g, false, true));
	// true indicates revese complement (second read)
	EXPECT_EQ(getStartKmerPos(k, testRead, g),
		getStartKmerPos(k, testRead, g, true));
}

TEST_F(GetStartKmerPosTest, EqualLengthMatchRegions)
{
	const string& seq = testRead.seq;
	BloomFilter bloom(bloomFilterSize);
	DBGBloom g(bloom);

	// This loop creates kmer match vector 011011
	for (unsigned i = 0; i < seq.length() - k + 1; i++) {
		// non-matching kmers
		if (i == 0 || i == 3)
			continue;
		bloom.loadSeq(k, seq.substr(i,k));
	}

	EXPECT_EQ(1U, getStartKmerPos(k, testRead, g, false, true));
	// true indicates revese complement (second read)
	EXPECT_EQ(1U, getStartKmerPos(k, testRead, g, true, true));
}

class CorrectSingleBaseErrorTest : public testing::Test {

protected:

	FastaRecord correctRead;
	FastaRecord singleErrorRead;
	string simulatedFalsePositive;
	static const int k = 6;
	static const size_t errorPos1 = 4;
	// set this large to avoid false positives
	static const int bloomFilterSize = 1000000;

	CorrectSingleBaseErrorTest() {
		Kmer::setLength(k);
		correctRead.seq = "TACAGTGCC";
		singleErrorRead.seq = correctRead.seq;
		singleErrorRead.seq[errorPos1] = 'C';
		simulatedFalsePositive = "TGCAGT";
	}

};


TEST_F(CorrectSingleBaseErrorTest, SingleError)
{
	BloomFilter bloom(bloomFilterSize);
	bloom.loadSeq(k, correctRead.seq);
	DBGBloom g(bloom);

	size_t correctedPos = numeric_limits<size_t>::max();
	FastaRecord read = singleErrorRead;
	bool success = correctSingleBaseError(g, k, read, correctedPos);

	ASSERT_TRUE(success);
	EXPECT_EQ(read.seq, read.seq);
	EXPECT_TRUE(correctedPos == errorPos1);
}

TEST_F(CorrectSingleBaseErrorTest, ReverseComplement)
{
	BloomFilter bloom(bloomFilterSize);
	bloom.loadSeq(k, correctRead.seq);
	DBGBloom g(bloom);

	size_t correctedPos = numeric_limits<size_t>::max();
	FastaRecord read = singleErrorRead;
	bool success = correctSingleBaseError(g, k, read, correctedPos, true);

	ASSERT_TRUE(success);
	EXPECT_EQ(read.seq, read.seq);
	EXPECT_TRUE(correctedPos == errorPos1);
}

TEST_F(CorrectSingleBaseErrorTest, NoError)
{
	BloomFilter bloom(bloomFilterSize);
	bloom.loadSeq(k, singleErrorRead.seq);
	DBGBloom g(bloom);

	size_t correctedPos = numeric_limits<size_t>::max();
	FastaRecord read = singleErrorRead;
	bool success = correctSingleBaseError(g, k, read, correctedPos);
	ASSERT_FALSE(success);
}

TEST_F(CorrectSingleBaseErrorTest, SkipFalsePositive)
{
	BloomFilter bloom(bloomFilterSize);
	bloom.loadSeq(k, correctRead.seq);
	bloom.loadSeq(k, simulatedFalsePositive);
	DBGBloom g(bloom);

	size_t correctedPos = numeric_limits<size_t>::max();
	FastaRecord read = singleErrorRead;
	bool success = correctSingleBaseError(g, k, read, correctedPos);
	ASSERT_TRUE(success);
	EXPECT_EQ(read.seq, read.seq);
	EXPECT_TRUE(correctedPos == errorPos1);
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
