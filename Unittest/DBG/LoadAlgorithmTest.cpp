#include "Assembly/SequenceCollection.h"
#include "Assembly/DBG.h"
#include "Assembly/AssemblyAlgorithms.h"
#include "Assembly/Options.h"
#include "Common/UnorderedSet.h"

#include <gtest/gtest.h>
#include <string>
#include <iostream>

using namespace std;

TEST(LoadAlgorithmTest, base)
{
	typedef SequenceCollectionHash Graph;
	Graph g;

	opt::kmerSize = 5;
	Kmer::setLength(5);

	Sequence seq("TAATGCCA");

	AssemblyAlgorithms::loadSequence(&g, seq);

	unordered_set<Kmer> kmers, expectedKmers;

	expectedKmers.insert(Kmer("TAATG"));
	expectedKmers.insert(Kmer("AATGC"));
	expectedKmers.insert(Kmer("ATGCC"));
	expectedKmers.insert(Kmer("TGCCA"));

	for (Graph::const_iterator it = g.begin();
			it != g.end(); ++it) {
		Kmer kmer(it->first);
#if 0
cerr << "visiting Kmer: " << kmer << "\n";
#endif
		ASSERT_TRUE(expectedKmers.find(kmer) != expectedKmers.end());
		expectedKmers.erase(kmer);
	}

	ASSERT_TRUE(expectedKmers.empty());
}
