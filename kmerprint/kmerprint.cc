/**
 * Print a kmer file. A kmer file is a serialized Google sparsehash.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */

#include "SequenceCollectionHash.h"
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
	if (opt::sequence)
		cout << seq.decode() << '\t';
	if (opt::adj)
		cout << seq.getExtension(SENSE) << '\t'
			<< seq.getExtension(ANTISENSE) << '\t';
	cout << seq.getMultiplicity() << '\n';
}

static void print(const PackedSeq& seq, extDirection sense)
{
	if (opt::sequence) {
		if (sense) {
			PackedSeq rc(seq);
			rc.reverseComplement();
			cout << rc.decode();
		} else
			cout << seq.decode();
		cout << '\t';
	}
	if (opt::adj)
		cout << seq.getExtension(sense) << '\t';
	cout << seq.getMultiplicity(sense) << '\n';
}

static void printFile(const char* path)
{
	SequenceCollectionHash c;
	c.load(path);
	for (SequenceCollectionHash::const_iterator it = c.begin();
			it != c.end(); ++it) {
		if (it->deleted())
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
