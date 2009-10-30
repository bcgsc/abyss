/**
 * Print a kmer file. A kmer file is a serialized Google sparsehash.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */

#include "SequenceCollection.h"
#include <algorithm>
#include <cassert>
#include <iostream>

using namespace std;

namespace opt {
	bool strands;
	bool sequence;
	bool adj;
}

static void print(const PackedSeq& seq)
{
	const KmerData& data = seq.second;
	if (opt::sequence)
		cout << seq.first.decode() << '\t';
	if (opt::adj)
		cout << data.getExtension(SENSE) << '\t'
			<< data.getExtension(ANTISENSE) << '\t';
	cout << data.getMultiplicity() << '\n';
}

static void print(const PackedSeq& seq, extDirection sense)
{
	const KmerData& data = seq.second;
	if (opt::sequence) {
		if (sense)
			cout << reverseComplement(seq.first).decode();
		else
			cout << seq.first.decode();
		cout << '\t';
	}
	if (opt::adj)
		cout << data.getExtension(sense) << '\t';
	cout << data.getMultiplicity(sense) << '\n';
}

static void printFile(const char* path)
{
	SequenceCollectionHash c;
	c.load(path);
	for (SequenceCollectionHash::const_iterator it = c.begin();
			it != c.end(); ++it) {
		if (it->second.deleted())
			continue;
		if (opt::strands) {
			print(*it, SENSE);
			print(*it, ANTISENSE);
		} else
			print(*it);
	}
}

int main(int argc, const char* argv[])
{
	assert(argc > 1);
	for_each(argv + 1, argv + argc, printFile);
	return 0;
}
