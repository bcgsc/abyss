#include "AssemblyAlgorithms.h"
#include "CommonUtils.h"
#include "DotWriter.h"
#include "ISequenceCollection.h"
#include "Options.h"
#include "SequenceCollectionHash.h"
#include "Timer.h"
#include <algorithm>
#include <cstdio>
#include <fstream>

static void popBubbles(/*const*/ ISequenceCollection* pSC)
{
	if (opt::bubbles <= 0)
		return;
	puts("Popping bubbles");
	unsigned totalPopped = 0;
	int i;
	for (i = 0; i < opt::bubbles; i++) {
		unsigned numPopped = AssemblyAlgorithms::popBubbles(pSC,
				opt::kmerSize);
		if (numPopped == 0)
			break;
		totalPopped += numPopped;
	}
	printf("Removed %d bubbles in %d rounds\n",
			totalPopped, i);
}

static void write_graph(const std::string& path,
		/*const*/ ISequenceCollection& c)
{
	if (path.length() == 0)
		return;
	printf("Writing graph to %s\n", path.c_str());
	ofstream out(path.c_str());
	DotWriter::write(out, c);
}

int main(int argc, char* const* argv)
{	
	Timer timer("Abyss");

	// Set stdout to be line buffered.
	setvbuf(stdout, NULL, _IOLBF, 0);

	opt::parse(argc, argv);
	
	// Load the phase space
	SequenceCollectionHash* pSC = new SequenceCollectionHash();
	//SequenceCollection* pSC = new SequenceCollection();
	
	for_each(opt::inFiles.begin(), opt::inFiles.end(),
			bind1st(ptr_fun(AssemblyAlgorithms::loadSequences), pSC));

	printf("total sequences: %d\n", pSC->count());

	printf("finalizing\n");
	pSC->finalize();

	AssemblyAlgorithms::generateAdjacency(pSC);
	

	if (opt::erode > 0) {
		puts("Eroding");
		unsigned totalEroded = 0;
		for (int i = 0; i < opt::erode; i++)
			totalEroded += AssemblyAlgorithms::erodeEnds(pSC);
		printf("Eroded %d tips in total\n",
				totalEroded);
	}

	int startTrimLength = 2;
	AssemblyAlgorithms::performTrim(pSC, startTrimLength);
		
	popBubbles(pSC);
	
	// Perform an additional trim at the max trim length to get rid of any new dead ends that formed during the bubble popping
	// These dead ends can happen when there are two overlapping bubbles and the second one is trimmed first (the bubble with only 2 branches)
	// There may be a better way to deal with this situation but this will suffice for the moment
	if (opt::bubbles > 0 && opt::trimLen > 0)
		while(AssemblyAlgorithms::trimSequences(pSC, opt::trimLen));

	write_graph(opt::graphPath, *pSC);

	AssemblyAlgorithms::splitAmbiguous(pSC);
	
	FastaWriter writer(opt::contigsPath.c_str());

	AssemblyAlgorithms::assemble(pSC, opt::readLen, opt::kmerSize, &writer);


	delete pSC;
	
	return 0;
}
