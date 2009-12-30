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
#include <cstdio>
#include <fstream>

using namespace std;

static void splitAmbiguousEdges(ISequenceCollection* pSC)
{
	unsigned marked = AssemblyAlgorithms::markAmbiguous(pSC);
	unsigned split = AssemblyAlgorithms::splitAmbiguous(pSC);
	assert(marked == split);
	(void)marked;
	(void)split;
}

static void removeLowCoverageContigs(ISequenceCollection* pSC)
{
	splitAmbiguousEdges(pSC);

	printf("Removing low-coverage contigs "
			"(mean k-mer coverage < %f)\n", opt::coverage);

	AssemblyAlgorithms::assemble(pSC);

	pSC->wipeFlag(SeqFlag(SF_MARK_SENSE | SF_MARK_ANTISENSE));
	opt::coverage = 0;
}

static void popBubbles(ISequenceCollection* pSC)
{
	puts("Popping bubbles");
	ofstream out;
	AssemblyAlgorithms::openBubbleFile(out);
	unsigned numPopped = AssemblyAlgorithms::popBubbles(pSC, out);
	assert(out.good());
	printf("Removed %d bubbles\n", numPopped);
}

static void write_graph(const string& path,
		const ISequenceCollection& c)
{
	if (path.empty())
		return;
	printf("Writing graph to %s\n", path.c_str());
	ofstream out(path.c_str());
	DotWriter::write(out, c);
}

int main(int argc, char* const* argv)
{
	Timer timer("Total");

	// Set stdout to be line buffered.
	setvbuf(stdout, NULL, _IOLBF, 0);

	opt::parse(argc, argv);

	SequenceCollectionHash* pSC = new SequenceCollectionHash();

	for_each(opt::inFiles.begin(), opt::inFiles.end(),
			bind1st(ptr_fun(AssemblyAlgorithms::loadSequences), pSC));
	printf("Loaded %zu k-mer\n", pSC->count());
	pSC->printLoad();
	assert(pSC->count() > 0);

	AssemblyAlgorithms::determineMinimumCoverage(
			AssemblyAlgorithms::coverageHistogram(*pSC));

generate_adjacency:
	puts("Generating adjacency");
	AssemblyAlgorithms::generateAdjacency(pSC);

	if (opt::erode > 0) {
		puts("Eroding tips");
		AssemblyAlgorithms::erodeEnds(pSC);
		assert(AssemblyAlgorithms::erodeEnds(pSC) == 0);
		pSC->cleanup();
		pSC->printLoad();
	}

	AssemblyAlgorithms::performTrim(pSC);

	if (opt::coverage > 0) {
		removeLowCoverageContigs(pSC);
		goto generate_adjacency;
	}

	if (opt::bubbles > 0)
		popBubbles(pSC);

	write_graph(opt::graphPath, *pSC);

	splitAmbiguousEdges(pSC);

	FastaWriter writer(opt::contigsPath.c_str());

	unsigned nContigs = AssemblyAlgorithms::assemble(pSC, &writer);
	if (nContigs == 0) {
		fputs("error: no contigs assembled\n", stderr);
		exit(EXIT_FAILURE);
	}

	delete pSC;

	return 0;
}
