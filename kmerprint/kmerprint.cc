/**
 * Print a kmer file. A kmer file is a serialized Google sparsehash.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */

#include "SequenceCollectionHash.h"
#include <cassert>
#include <iostream>

using namespace std;

namespace opt {
	bool sequence = true;
	bool strands = true;
}

static void print(const PackedSeq& seq)
{
	if (opt::sequence)
		cout << seq.decode() << '\t';
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
	cout << seq.getMultiplicity(sense) << '\n';
}

int main(int argc, const char* argv[])
{
	SequenceCollectionHash c;
	assert(argc > 1);
	c.load(argv[1]);
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
	return 0;
}
