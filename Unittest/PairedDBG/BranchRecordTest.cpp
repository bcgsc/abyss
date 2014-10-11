#include "PairedDBG/BranchRecord.h"
#include "Common/Options.h"
#include "PairedDBG/Options.h"

#include <gtest/gtest.h>
#include <vector>

using namespace std;

TEST(BranchRecordTest, Sequence)
{
	// length of each kmer in kmer pair
	Kmer::setLength(2);
	// space between kmer pair
	unsigned delta = 2;
	// the length of both kmers plus the gap
	KmerPair::setLength(Kmer::length() * 2 + delta);

	// sequence for branch: TAGGGATT	
	// kmer pairs:          TA  GA
    //                       AG  AT
	//                        GG  TT

	std::pair<KmerPair,KmerPairData>
		kmerPair1(KmerPair("TA", "GA"), KmerPairData()),
		kmerPair2(KmerPair("AG", "AT"), KmerPairData()),
		kmerPair3(KmerPair("GG", "TT"), KmerPairData());
		
	// test sequence reconstruction in forward dir

	BranchRecord forwardBranch(SENSE);
	forwardBranch.push_back(kmerPair1);
	forwardBranch.push_back(kmerPair2);
	forwardBranch.push_back(kmerPair3);
	ASSERT_EQ("TAGGGATT", (Sequence)forwardBranch);

	// test sequence reconstruction in reverse dir

	BranchRecord reverseBranch(ANTISENSE);
	reverseBranch.push_back(kmerPair3);
	reverseBranch.push_back(kmerPair2);
	reverseBranch.push_back(kmerPair1);
	ASSERT_EQ("TAGGGATT", (Sequence)reverseBranch);
	
	// test sequence reconstruction with "N"s (forward dir)

	BranchRecord shortBranchForward(SENSE);
	shortBranchForward.push_back(kmerPair1);
	shortBranchForward.push_back(kmerPair2);
	ASSERT_EQ("TAGNGAT", (Sequence)shortBranchForward);
	
	// test sequence reconstruction with "N"s (reverse dir)

	BranchRecord shortBranchReverse(ANTISENSE);
	shortBranchReverse.push_back(kmerPair2);
	shortBranchReverse.push_back(kmerPair1);
	ASSERT_EQ("TAGNGAT", (Sequence)shortBranchReverse);
}
