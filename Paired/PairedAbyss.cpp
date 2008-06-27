#include <stdio.h>
#include <math.h>
#include <iostream>
#include "PairedAbyss.h"
#include "SequenceCollectionHash.h"
#include "AssemblyAlgorithms.h"
#include "PairedAlgorithms.h"
#include "Options.h"
#include "FastaReader.h"
#include "VisitAlgorithms.h"
#include "Stats.h"

int main(int argc, char** argv)
{
	// Read in a trimmed.fa file of reads, a contigs file (from single-end abyss), the file with the original reads in it and the assembly parameters
	std::string finalReadsFile(argv[1]);
	std::string readsFile(argv[2]);
	std::string contigFile(argv[3]);
	int readLen = atoi(argv[4]);
	int kmer = atoi(argv[5]);

	
	// kinda hacky
	opt::kmerSize = kmer;
	opt::readLen = readLen;
	
	// Load the adjacency information from the trimmed sequence file
	SequenceCollectionHash* pSC = new SequenceCollectionHash();
	AssemblyAlgorithms::loadSequences(pSC, finalReadsFile);
	pSC->finalize();
	AssemblyAlgorithms::generateAdjacency(pSC);	
	
	// Read in the initial contigs
	ContigMap contigMap;
	PairedAlgorithms::readContigMap(contigFile, contigMap);
	
	// Generate the initial alignment DB, it will be populated during the generation of the graph
	AlignmentCache alignDB;
	
	// Generate the contig graph from the adjacency information and the input contigs
	ContigGraph* pContigGraph = new ContigGraph();
	PairedAlgorithms::generateGraph(pContigGraph, contigMap, pSC, kmer, &alignDB);
	
	// Add the pairs of each sequence onto the contigs they belong to

	// Read in the pairings
	PairRecord pairRecord;
	
	if(argc <= 6)
	{
		LoadPairsRecord(readsFile, kmer, pairRecord);
		pairRecord.serialize("PairsCache.bin");
	}
	else
	{
		pairRecord.unserialize(argv[6]);	
	}
		
	// Validate the db
	pContigGraph->iterativeVisit(DBValidator(&alignDB));	
	printf("Database validated\n");
	
	printf("Linking pairs to contigs\n");
	// Add the pairs to contigs
	pContigGraph->iterativeVisit(PairAdder(&pairRecord));
	
	printf("Setting up probability models\n");
	// Generate the histogram of pair distances
	Histogram pairDistHist;
	pContigGraph->iterativeVisit(PairedDistHistogramGenerator(&pairDistHist, 300));
	
	PairedStats pairedStats;
	pairedStats.generateStats(pairDistHist);
	
	// Generate the paired resolver policy based on the empirical pdf
	PairedResolvePolicy pairedResolvePolicy(pairedStats.getPDF());
	
	printf("Resolving initial self pairings\n");
	// Resolve all the self-pairs
	pContigGraph->iterativeVisit(PairedResolveVisitor(&pairedResolvePolicy));
	
	PairedMerger pairedMerger(kmer, pairedStats.getMax(), &pairedResolvePolicy, &alignDB);
	/*
	pContigGraph->printVertex("507");
	pContigGraph->printVertex("11");
	pContigGraph->printVertex("805");
	pContigGraph->printVertex("2419");
	pContigGraph->printVertex("1747");
	pContigGraph->printVertex("1750");
	pContigGraph->printVertex("1008");
	pContigGraph->printVertex("910");
	pContigGraph->printVertex("936");
	pContigGraph->printVertex("1051");
	return 1;
	*/

	
	ContigID id = "2085";
	//ContigID id = "2033";
	//ContigID id = "2634";
	SequenceDataCost dataCost(kmer);
	//while(pairedMerger.resolve(pContigGraph, id))
	for(int i = 0; i < 1; ++i)
	{
		pairedMerger.resolve(pContigGraph, id);
	}
	pContigGraph->getDataForVertex(id).printPairAlignments(ANTISENSE, PSF_ALL);
	pContigGraph->printVertex(id);
	return 1;
	
	//pContigGraph->reducePaired(pairedMerger);
	//data.printPairAlignments(SENSE);
	//const_cast<ContigData*>(&data)->computeTestStat(pairedStats.GetPDF());
	
	FastaWriter* pWriter = new FastaWriter("PostPaired.fa");
	pContigGraph->iterativeVisit(ContigDataOutputter(pWriter));
	delete pWriter;
	
	delete pSC;
	delete pContigGraph;
	
	printf("Done\n");
}

void LoadPairsRecord(std::string file, int kmerSize, PairRecord& pairRecord)
{
	// Read in the fasta file
	FastaReader reader(file.c_str());
	int count = 0;
	while(reader.isGood())
	{
		PackedSeq seq1 = reader.ReadSequence();
		PackedSeq seq2 = reader.ReadSequence();
				
		int len1 = seq1.getSequenceLength();
		int len2 = seq2.getSequenceLength();
		
		assert(kmerSize <= len1);
		assert(len1 == len2);
		
		for(int i = 0; i < len1 - kmerSize  + 1; i++)
		{
			PackedSeq sub1 = seq1.subseq(i, kmerSize);
			
			// Read the pair from the end inwards
			PackedSeq sub2 = seq2.subseq((len2 - i) - kmerSize, kmerSize);

			pairRecord.addPairs(sub1, sub2);
			count++;
		}
	}
	
	printf("Found %d kmer pairs\n", count);
}
