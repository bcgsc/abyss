#include <stdio.h>

#include <vector>
#include "Sequence.h"
#include "Reader.h"
#include "PathDriver.h"
#include "PairRecord.h"
#include "Writer.h"
#include "SeqRecord.h"

int main(int argv, char** argc)
{
	Reader fileReader;
	
	if(argv < 2)
	{
		printf("usage: PathWalker <fasta file> <insert length>\n");
		exit(1);
	}
	
	const char* fastaFile = argc[1];
	int insertLength = atoi(argc[2]);
#if 0
	PSequenceVector allSeqs;
	bool result = fileReader.readFasta(fastaFile, allSeqs);
	
	if(result != true && allSeqs.size() > 0)
	{
		printf("error reading fasta file, exiting\n");
	}
	
	// Create the pair mapping from all the seqs
	// pairs are assumed to be in contiguous entries in the array
	PairRecord pairRecord(allSeqs);

	// Get the length of the reads
	// The read length is assumed to be uniform for all input sequences
	int readLength = allSeqs.front().length();
	
	//printf("read %d sequences of length %d\n", allSeqs.size(), readLength);	
	
	// Create phase space
	PhaseSpace phaseSpace(readLength);

	phaseSpace.addReads(allSeqs);
	
	// create file writer
	Writer writer("contigs.fa");
		
	// Create the master record of sequences
	SeqRecord multiplicityRecord;
	multiplicityRecord.addMultipleSequences(allSeqs);
	
	SeqRecord extensionRecord;
	
	int totalSeqs = allSeqs.size() - 1;
	
	for(int i = 0; i < allSeqs.size(); i++)
	{
		Sequence seedSeq = allSeqs[i];

		if(extensionRecord.getMultiplicity(seedSeq) > 0)
		{
			continue;	
		}	
	
		PathDriver driver(seedSeq, SENSE, &phaseSpace, &pairRecord, &multiplicityRecord, &extensionRecord);
		Path finalPathSense = driver.run();
			
		PathDriver driver2(seedSeq, ANTISENSE, &phaseSpace, &pairRecord, &multiplicityRecord, &extensionRecord);
		Path finalPathAntisense = driver2.run();
		
		finalPathAntisense.mergePath(finalPathSense, false, true, true);
		
		finalPathAntisense.print();
		
		if(finalPathAntisense.getSequence().length() >= 100)
		{
			writer.writeContig(finalPathAntisense.getSequence().c_str());
		}
	}
		
	
	/*
	for(int index = 0; index < 1500000; index++)
	{
		// this is weak 
		int randNum = rand();
		int randomIndex = ((double)randNum / (double)RAND_MAX) * (double)totalSeqs;
		printf("rand %d randidx %d randmax %d\n", randNum, randomIndex, RAND_MAX);
		Sequence seedSeq = allSeqs[randomIndex];
		
		if(extensionRecord.getMultiplicity(seedSeq) > 0)
		{
			continue;	
		}
		

		PathDriver driver(seedSeq, SENSE, &phaseSpace, &pairRecord, &multiplicityRecord, &extensionRecord);
		Path finalPathSense = driver.run();
		
		PathDriver driver2(seedSeq, ANTISENSE, &phaseSpace, &pairRecord, &multiplicityRecord, &extensionRecord);
		Path finalPathAntisense = driver2.run();
		
		finalPathAntisense.mergePath(finalPathSense, false, true, true);
		
		if(finalPathAntisense.getSequence().length() > 200)
		{
			writer.writeContig(finalPathAntisense.getSequence().c_str());
		}
	}
	*/
#endif	
	return 0;
}
