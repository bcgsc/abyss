#include "Konnector/DBGBloomAlgorithms.h"
#include "Bloom/Bloom.h"
#include "Bloom/CascadingBloomFilter.h"
#include "Common/Sequence.h"
#include <gtest/gtest.h>
#include <string>

using namespace std;

/*
 * Tests for getStartKmerPos() function, which does
 * the following:
 *
 * Choose a suitable starting kmer for a path search and
 * return its position. More specifically, find the kmer
 * closest to the end of the given sequence that is followed by
 * at least (numMatchesThreshold - 1) consecutive kmers that
 * are also present in the Bloom filter de Bruijn graph. If there
 * is no sequence of matches of length numMatchesThreshold,
 * use the longest sequence of matching kmers instead.
 *
 * @param seq sequence in which to find start kmer
 * @param k kmer size
 * @param g de Bruijn graph
 * @param numMatchesThreshold if we encounter a sequence
 * of numMatchesThreshold consecutive kmers in the Bloom filter,
 * choose the kmer at the beginning of that sequence
 * @return position of chosen start kmer
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
	Bloom::loadSeq(bloom, k, testRead.seq);
	DBGBloom<BloomFilter> g(bloom);
	const unsigned numMatchesThreshold = 1;

	EXPECT_EQ(5U, getStartKmerPos(testRead, k, FORWARD, g,
		numMatchesThreshold));
}

TEST_F(GetStartKmerPosTest, FullReadMismatch)
{
	BloomFilter bloom(bloomFilterSize);
	// Leave the bloom filter empty to generate a mismatch
	// for every kmer in the read.
	DBGBloom<BloomFilter> g(bloom);
	EXPECT_EQ(NO_MATCH, getStartKmerPos(testRead, k, FORWARD, g));
}

TEST_F(GetStartKmerPosTest, NumMatchesThreshold)
{
	const string& seq = testRead.seq;
	BloomFilter bloom(bloomFilterSize);
	DBGBloom<BloomFilter> g(bloom);

	// This loop creates kmer match vector 101101
	for (unsigned i = 0; i < seq.length() - k + 1; i++) {
		// non-matching kmers
		if (i == 1 || i == 4)
			continue;
		Bloom::loadSeq(bloom, k, seq.substr(i,k));
	}

	unsigned numMatchesThreshold;

	numMatchesThreshold = 1;
	EXPECT_EQ(5U, getStartKmerPos(testRead, k, FORWARD, g,
		numMatchesThreshold));

	numMatchesThreshold = 2;
	EXPECT_EQ(2U, getStartKmerPos(testRead, k, FORWARD, g,
		numMatchesThreshold));

	numMatchesThreshold = 3;
	EXPECT_EQ(2U, getStartKmerPos(testRead, k, FORWARD, g,
		numMatchesThreshold));

	numMatchesThreshold = 1;
	EXPECT_EQ(0U, getStartKmerPos(testRead, k, REVERSE, g,
		numMatchesThreshold));

	numMatchesThreshold = 2;
	EXPECT_EQ(3U, getStartKmerPos(testRead, k, REVERSE, g,
		numMatchesThreshold));

	numMatchesThreshold = 3;
	EXPECT_EQ(3U, getStartKmerPos(testRead, k, REVERSE, g,
		numMatchesThreshold));
}

TEST_F(GetStartKmerPosTest, EqualLengthMatchRegions)
{
	const string& seq = testRead.seq;
	BloomFilter bloom(bloomFilterSize);
	DBGBloom<BloomFilter> g(bloom);

	// This loop creates kmer match vector 011011
	for (unsigned i = 0; i < seq.length() - k + 1; i++) {
		// non-matching kmers
		if (i == 0 || i == 3)
			continue;
		Bloom::loadSeq(bloom, k, seq.substr(i,k));
	}

	unsigned numMatchesThreshold;

	numMatchesThreshold = 2;
	EXPECT_EQ(4U, getStartKmerPos(testRead, k, FORWARD, g,
		numMatchesThreshold));

	numMatchesThreshold = 2;
	EXPECT_EQ(2U, getStartKmerPos(testRead, k, REVERSE, g,
		numMatchesThreshold));
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
	Bloom::loadSeq(bloom, k, correctRead.seq);
	DBGBloom<BloomFilter> g(bloom);

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
	Bloom::loadSeq(bloom, k, correctRead.seq);
	DBGBloom<BloomFilter> g(bloom);

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
	Bloom::loadSeq(bloom, k, singleErrorRead.seq);
	DBGBloom<BloomFilter> g(bloom);

	size_t correctedPos = numeric_limits<size_t>::max();
	FastaRecord read = singleErrorRead;
	bool success = correctSingleBaseError(g, k, read, correctedPos);
	ASSERT_FALSE(success);
}

TEST_F(CorrectSingleBaseErrorTest, SkipFalsePositive)
{
	BloomFilter bloom(bloomFilterSize);
	Bloom::loadSeq(bloom, k, correctRead.seq);
	Bloom::loadSeq(bloom, k, simulatedFalsePositive);
	DBGBloom<BloomFilter> g(bloom);

	size_t correctedPos = numeric_limits<size_t>::max();
	FastaRecord read = singleErrorRead;
	bool success = correctSingleBaseError(g, k, read, correctedPos);
	ASSERT_TRUE(success);
	EXPECT_EQ(read.seq, read.seq);
	EXPECT_TRUE(correctedPos == errorPos1);
}
