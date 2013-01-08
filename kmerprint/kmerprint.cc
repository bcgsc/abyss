/**
 * Print a kmer file. A kmer file is a serialized Google sparsehash.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */

#include "SequenceCollection.h"
#include "Uncompress.h"
#include <algorithm>
#include <cassert>
#include <getopt.h>
#include <iostream>
#include <sstream>

using namespace std;

#define PROGRAM "kmerprint"

namespace opt {
	static unsigned k;
	static bool strands;
}

static const struct option longopts[] = {
	{ "kmer",    required_argument, NULL, 'k' },
};

static const char shortopts[] = "k:";

static void print(const SequenceCollectionHash::value_type& seq)
{
	const KmerData& data = seq.second;
	cout << seq.first.str()
		<< '\t' << data.getMultiplicity()
		<< '\t' << data.getExtension(SENSE)
		<< '\t' << data.getExtension(ANTISENSE)
		<< '\n';
}

static void print(const SequenceCollectionHash::value_type& seq,
		extDirection sense)
{
	const KmerData& data = seq.second;
	cout << (sense ? reverseComplement(seq.first).str()
			: seq.first.str())
		<< '\t' << data.getMultiplicity(sense)
		<< '\t' << data.getExtension(sense)
		<< '\n';
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

int main(int argc, char* argv[])
{
	assert(argc > 1);

	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': exit(EXIT_FAILURE);
			case 'k': arg >> opt::k; break;
		}
	}

	if (opt::k <= 0) {
		cerr << PROGRAM ": " << "missing -k,--kmer option\n";
		exit(EXIT_FAILURE);
	}

	Kmer::setLength(opt::k);

	for_each(argv + optind, argv + argc, printFile);
	return 0;
}
