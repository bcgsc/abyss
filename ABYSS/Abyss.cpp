#include <stdio.h>

#include <vector>
#include <stdio.h>
#include <deque>
#include <iostream>
#include <fstream>
#include "Abyss.h"
#include "CommonUtils.h"
#include "DotWriter.h"
#include "SequenceCollection.h"
#include "SequenceCollectionHash.h"
#include "AssemblyAlgorithms.h"

int main(int argc, char** argv)
{	

	if(argc < 4 || argv[1] == "--help")
	{
		printUsage();
		exit(1);
	}

	std::string fastaFile = argv[1];
	int readLen = atoi(argv[2]);
	int kmerSize = atoi(argv[3]);
	
	bool noTrim = false;
	if(argc == 5)
	{
		std::string flags = argv[4];
		if(flags == "-notrim")
		{
			noTrim = true;
		}
	}

	// Load the phase space
	SequenceCollectionHash* pSC = new SequenceCollectionHash();
	//SequenceCollection* pSC = new SequenceCollection();
	
	loadSequences(pSC, fastaFile, readLen, kmerSize);

	printf("total sequences: %d\n", pSC->count());
	
	printf("finalizing\n");
	pSC->finalize();

	generateAdjacency(pSC);

	if(!noTrim)
	{
		performTrim(pSC, readLen, kmerSize);
	}
	
	outputSequences("trimmed.fa", pSC);
	
	// Remove bubbles
	popBubbles(pSC, kmerSize);

	puts("Building graph.dot...");
	ofstream dot_out("graph.dot");
	DotWriter::write(dot_out, *pSC);
	puts("Done.");

	splitAmbiguous(pSC);
	
	assemble(pSC, readLen, kmerSize);

	delete pSC;
	
	return 0;
}

/*
void outputBranchSizes(AssemblyData* pSS, Coord4 minCoord, Coord4 maxCoord)
{	 
	static int trimNum = 0;
	static std::ofstream branchLog("branchLog.txt");
	//printf("seqs before trimming: %d\n", pSS->countAll());
	int numBranchesRemoved = 0;
	int count = 0;
	for(SequenceCollectionIter iter = pSS->getStartIter(); iter != pSS->getEndIter(); iter++)
	{	
		if(iter->isFlagSet(SF_DELETE))
		{
			continue;
		}
				
		// does this sequence have an parent/child extension?
		bool hasChild = pSS->hasChild(*iter);
		bool hasParent = pSS->hasParent(*iter);
		extDirection dir;
		
		if(!hasChild && !hasParent)
		{
			// remove this sequence, it has no extensions
			pSS->removeSequence(*iter);
			continue;
		}
		else if(!hasChild)
		{
			dir = ANTISENSE;
		}
		else if(!hasParent)
		{
			dir = SENSE;
		}
		else
		{
			continue;
		}

		PSequenceVector branchElements;
		extDirection oppositeDir = oppositeDirection(dir);
					
		PackedSeq currSeq = *iter;
		SeqRecord loopCheck;			
		bool stop = false;
		
		while(!stop)
		{				
			HitRecord hr = pSS->calculateExtension(currSeq, dir);
			HitRecord oppHr = pSS->calculateExtension(currSeq, oppositeDir);
				
			if(oppHr.getNumHits() > 1)
			{
				//printf("stopped because of reverse branch\n");
				stop = true;	
			}
			else if(hr.getNumHits() == 0 || hr.getNumHits() > 1)
			{
				branchElements.push_back(currSeq);
				//printf("stopped because of noext/ambi branch\n");
				stop = true;
			}
			else
			{
				branchElements.push_back(currSeq);
					
				// good ext
				currSeq = hr.getFirstHit().seq;
			}
		}
		branchLog << trimNum << " " << branchElements.size() << std::endl;
	}
	
	trimNum++;
}

int mainHash(int argc, char** argv)
{
	std::string fastaFile = argv[1];
	int readLen = atoi(argv[2]);
	int kmerSize = atoi(argv[3]);
	FastaReader* reader = new FastaReader(fastaFile.c_str());
	
	int NUM_PROCS = 10000000;
	int* vals = new int[NUM_PROCS];
	for(int i = 0; i < NUM_PROCS; i++)
	{
		vals[i] = 0;
	}
	
	int count = 0;
	while(reader->isGood())
	{
		PackedSeq seq = reader->ReadSequence();
		for(int i = 0; i < seq.getSequenceLength() - kmerSize  + 1; i++)
		{
			PackedSeq sub = seq.subseq(i, kmerSize);		
		
			unsigned int code = sub.getHashCode();

		
			unsigned int id = code % NUM_PROCS;
			vals[id]++;
			count++;
		}
		
		//printf("%d\n", id);
	}
	
	int max = 0;
	int min = count;
	int sum = 0;
	for(int i = 0; i < NUM_PROCS; i++)
	{
		if(vals[i] > max)
		{
			max = vals[i];
		}
		if(vals[i] < min)
		{
			min = vals[i];
		}
		
		sum += vals[i];
		
		printf("%d\n", vals[i]);
	}
	
	//printf("max: %d min: %d ratio: %lf mean: %lf total: %d\n", max, min, (double)min/(double)max, (double)sum/(double)NUM_PROCS, count);
	
	delete [] vals;
	delete reader;
	reader = NULL;
}
*/


void printUsage()
{
	printf("usage: ABYSS <reads fasta file> <max read length> <kmer size>\n");	
}
