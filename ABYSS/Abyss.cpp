#include <stdio.h>

#include <vector>
#include <stdio.h>
#include <deque>
#include <iostream>
#include <fstream>
#include "Abyss.h"
#include "CommonUtils.h"
#include "DotWriter.h"
#include "ISequenceCollection.h"
#include "Options.h"
#include "SequenceCollection.h"
#include "SequenceCollectionHash.h"
#include "AssemblyAlgorithms.h"
#include "Timer.h"

static void write_graph(const std::string& path,
		/*const*/ ISequenceCollection& c)
{
	cout << "Building " << path << "..." << endl;
	ofstream out(path.c_str());
	DotWriter::write(out, c);
	cout << "Done." << endl;
}

int main(int argc, char* const* argv)
{	
	Timer timer("Abyss");
	opt::parse(argc, argv);
	
	// Load the phase space
	SequenceCollectionHash* pSC = new SequenceCollectionHash();
	//SequenceCollection* pSC = new SequenceCollection();
	
	AssemblyAlgorithms::loadSequences(pSC, opt::inFile, opt::readLen, opt::kmerSize);

	printf("total sequences: %d\n", pSC->count());

	printf("finalizing\n");
	pSC->finalize();

	AssemblyAlgorithms::generateAdjacency(pSC);

	AssemblyAlgorithms::performTrim(pSC);
	
	AssemblyAlgorithms::outputSequences("trimmed.fa", pSC);
	write_graph("trimmed.dot", *pSC);

	// Remove bubbles
	while(AssemblyAlgorithms::popBubbles(pSC, opt::kmerSize));
	
	// Perform an additional trim at the max trim length to get rid of any new dead ends that formed during the bubble popping
	// These dead ends can happen when there are two overlapping bubbles and the second one is trimmed first (the bubble with only 2 branches)
	// There may be a better way to deal with this situation but this will suffice for the moment
	if (opt::trimLen > 0)
	{
		while(AssemblyAlgorithms::trimSequences(pSC, opt::trimLen));
	}

	write_graph("graph.dot", *pSC);

	AssemblyAlgorithms::splitAmbiguous(pSC);
	
	FastaWriter writer("contigs.fa");

	AssemblyAlgorithms::assemble(pSC, opt::readLen, opt::kmerSize, &writer);


	delete pSC;
	
	return 0;
}
