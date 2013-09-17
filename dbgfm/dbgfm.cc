/**
 * de Bruin Graph assembler using an FM-index
 * Copyright 2013 Shaun Jackman
 */

#include "Assembly/Options.h"
#include "Common/Timer.h"
#include "Common/Uncompress.h"
#include "DataLayer/FastaWriter.h"
#include "dbgfm/DBGAlgorithms.h"
#include "dbgfm/dbgfm.h"

#include <cstdio> // for setvbuf
#include <iostream>

using namespace std;

static void assemble(const string& pathIn, const string& pathOut)
{
	Timer timer(__func__);
	(void)pathIn;
	assert(pathIn.empty());
	assert(opt::inFiles.size() == 1);
	DBGFM g(opt::kmerSize, opt::inFiles[0]);
	cout << "Loaded " << num_vertices(g) << " k-mer\n";

	cout << "Generating adjacency" << endl;
	dbg::generateAdjacency(g);

	cout << "Marking ambiguous vertices" << endl;
	dbg::markAmbiguous(g);

	cout << "Assembling contigs" << endl;
	FastaWriter writer(pathOut.c_str());
	unsigned nContigs = dbg::assemble(g, &writer);
	if (nContigs == 0) {
		cerr << "error: no contigs assembled\n";
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char* const* argv)
{
	Timer timer("Total");

	// Set stdout to be line buffered.
	setvbuf(stdout, NULL, _IOLBF, 0);

	opt::parse(argc, argv);

	bool krange = opt::kMin != opt::kMax;
	if (krange)
		cout << "Assembling k=" << opt::kMin << "-" << opt::kMax
				<< ":" << opt::kStep << endl;
	assert(!krange);

	for (unsigned k = opt::kMin; k <= opt::kMax; k += opt::kStep) {
		if (krange)
			cout << "Assembling k=" << k << endl;
		opt::kmerSize = k;
		Kmer::setLength(k);

		if (k > opt::kMin) {
			// Reset the assembly options to defaults.
			opt::erode = (unsigned)-1;
			opt::erodeStrand = (unsigned)-1;
			opt::coverage = -1;
			opt::trimLen = k;
			opt::bubbleLen = 3*k;
		}

		ostringstream k0, k1;
		if (k > opt::kMin)
			k0 << "contigs-k" << k - opt::kStep << ".fa";
		if (k < opt::kMax)
			k1 << "contigs-k" << k << ".fa";
		else
			k1 << opt::contigsPath.c_str();
		assemble(k0.str(), k1.str());
	}
	return 0;
}
