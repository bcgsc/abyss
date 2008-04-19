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
	opt::parse(argc, argv);
	
	// Load the phase space
	SequenceCollectionHash* pSC = new SequenceCollectionHash();
	//SequenceCollection* pSC = new SequenceCollection();
	
	loadSequences(pSC, opt::fastaFile, opt::readLen, opt::kmerSize);

	printf("total sequences: %d\n", pSC->count());
	
	printf("finalizing\n");
	pSC->finalize();

	generateAdjacency(pSC);

	performTrim(pSC);
	
	outputSequences("trimmed.fa", pSC);
	write_graph("trimmed.dot", *pSC);

	// Remove bubbles
	popBubbles(pSC, opt::kmerSize);

	write_graph("graph.dot", *pSC);

	splitAmbiguous(pSC);
	
	assemble(pSC, opt::readLen, opt::kmerSize);

	delete pSC;
	
	return 0;
}
