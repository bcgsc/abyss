#include "Assembly/Options.h"
#include "AssemblyAlgorithms.h"
#include "DotWriter.h"
#include "FastaWriter.h"
#include "Histogram.h"
#include "ISequenceCollection.h"
#include "SequenceCollection.h"
#include "Timer.h"
#include "Uncompress.h"
#include <algorithm>
#include <cstdio> // for setvbuf
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

static void removeLowCoverageContigs(ISequenceCollection* pSC)
{
	AssemblyAlgorithms::markAmbiguous(pSC);

	cout << "Removing low-coverage contigs "
			"(mean k-mer coverage < " << opt::coverage << ")\n";
	AssemblyAlgorithms::assemble(pSC);
	AssemblyAlgorithms::splitAmbiguous(pSC);

	opt::coverage = 0;
}

static void popBubbles(ISequenceCollection* pSC)
{
	cout << "Popping bubbles" << endl;
	ofstream out;
	AssemblyAlgorithms::openBubbleFile(out);
	unsigned numPopped = AssemblyAlgorithms::popBubbles(pSC, out);
	assert(out.good());
	cout << "Removed " << numPopped << " bubbles\n";
}

static void write_graph(const string& path,
		const ISequenceCollection& c)
{
	if (path.empty())
		return;
	cout << "Writing graph to `" << path << "'\n";
	ofstream out(path.c_str());
	DotWriter::write(out, c);
}

static void assemble(const string& pathIn, const string& pathOut)
{
	Timer timer(__func__);
	SequenceCollectionHash* pSC = new SequenceCollectionHash();

	if (!pathIn.empty())
		AssemblyAlgorithms::loadSequences(pSC, pathIn.c_str());
	for_each(opt::inFiles.begin(), opt::inFiles.end(),
			bind1st(ptr_fun(AssemblyAlgorithms::loadSequences), pSC));
	cout << "Loaded " << pSC->count() << " k-mer\n";
	pSC->shrink();
	assert(pSC->count() > 0);

	AssemblyAlgorithms::setCoverageParameters(
			AssemblyAlgorithms::coverageHistogram(*pSC));

	cout << "Generating adjacency" << endl;
	AssemblyAlgorithms::generateAdjacency(pSC);

erode:
	if (opt::erode > 0) {
		cout << "Eroding tips" << endl;
		AssemblyAlgorithms::erodeEnds(pSC);
		assert(AssemblyAlgorithms::erodeEnds(pSC) == 0);
		pSC->cleanup();
	}

	AssemblyAlgorithms::performTrim(pSC);
	pSC->cleanup();

	if (opt::coverage > 0) {
		removeLowCoverageContigs(pSC);
		pSC->wipeFlag(SeqFlag(SF_MARK_SENSE | SF_MARK_ANTISENSE));
		pSC->cleanup();
		goto erode;
	}

	if (opt::bubbleLen > 0)
		popBubbles(pSC);

	write_graph(opt::graphPath, *pSC);

	AssemblyAlgorithms::markAmbiguous(pSC);

	FastaWriter writer(pathOut.c_str());
	unsigned nContigs = AssemblyAlgorithms::assemble(pSC, &writer);
	if (nContigs == 0) {
		cerr << "error: no contigs assembled\n";
		exit(EXIT_FAILURE);
	}

	delete pSC;
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

	for (int k = opt::kMin; k <= opt::kMax; k += opt::kStep) {
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
