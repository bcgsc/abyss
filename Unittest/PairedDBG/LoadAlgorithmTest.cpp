#include "PairedDBG/AssemblyAlgorithms.h"
#include "PairedDBG/Options.h"
#include "Common/UnorderedSet.h"

#include <gtest/gtest.h>
#include <string>
#include <iostream>

using namespace std;

TEST(LoadAlgorithmTest, base)
{
	SequenceCollectionHash g;

	// length of each kmer in kmer pair
	Kmer::setLength(2);
	// space between kmer pair
	unsigned delta = 2;
	// the length of both kmers plus the gap
	KmerPair::setLength(Kmer::length() * 2 + delta);

	// see test.png for an image of the paired
	// de Bruijn graph for this sequence
	Sequence seq("TAATGCCATGGGATGTT");

	AssemblyAlgorithms::loadSequence(&g, seq);

	unordered_set<KmerPair> kmerPairs, expectedKmerPairs;

	expectedKmerPairs.insert(KmerPair("TAGC"));
	expectedKmerPairs.insert(KmerPair("AACC"));
	expectedKmerPairs.insert(KmerPair("ATCA"));
	expectedKmerPairs.insert(KmerPair("GCTG"));
	expectedKmerPairs.insert(KmerPair("CCGG"));
	expectedKmerPairs.insert(KmerPair("CAGG"));
	expectedKmerPairs.insert(KmerPair("ATGA"));
	expectedKmerPairs.insert(KmerPair("GGTG"));
	expectedKmerPairs.insert(KmerPair("GGGT"));
	expectedKmerPairs.insert(KmerPair("GATT"));

	for (ISequenceCollection::iterator it = g.begin();
			it != g.end(); it++) {
		KmerPair kmerPair(it->first);
#if 0
cerr << "visiting KmerPair: " << kmerPair << "\n";
#endif
		ASSERT_TRUE(expectedKmerPairs.find(kmerPair) != expectedKmerPairs.end());
		expectedKmerPairs.erase(kmerPair);
	}

	ASSERT_TRUE(expectedKmerPairs.empty());
}
