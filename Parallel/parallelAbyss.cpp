#include <stdio.h>

#include <vector>
#include <stdio.h>
#include <mpi.h>
#include "parallelAbyss.h"
#include "CommonUtils.h"
#include "FastaReader.h"
#include "DistributedPhaseSpace.h"

int rank;

int main(int argc, char** argv)
{	
	
	if(argc < 4 || argv[1] == "--help")
	{
		printUsage();
		exit(1);
	}
	
	std::string fastaFile = argv[1];
	int kmerSize = atoi(argv[2]);
	//int numTrims = atoi(argv[3]);
	
	
   double starttime, endtime; 
   starttime = MPI::Wtime();
	
	// start mpi process
	MPI::Init(argc,argv);
	
	// get my rank
	rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	
	if(rank == 0)
	{
		printf("%d: starting read\n", rank);
		doRead(kmerSize, fastaFile, size - rank);
	}
	else
	{
		DistributedPhaseSpace dps(rank, kmerSize);
		dps.MessageLoop();
	}
	
	MPI::Finalize();
	
	endtime = MPI::Wtime();
	
	printf("%d elapsed seconds: %f\n", rank, endtime - starttime);
}

void doRead(int kmerSize, std::string file, int numDataNodes)
{
	FastaReader* reader = new FastaReader(file.c_str());
	CommLayer commLayer(0, kmerSize);
	
	// Load phase space
	int count = 0;

	PSequenceVector seqs;
	while(reader->isGood())
	{
		PackedSeq seq = reader->ReadSequence();		
		assert(kmerSize <= seq.getSequenceLength());
		
		for(int i = 0; i < seq.getSequenceLength() - kmerSize  + 1; i++)
		{
			PackedSeq sub = seq.subseq(i, kmerSize);
			
			// Send the sequence to the node
			commLayer.SendSequence(1, sub, APM_SEQLOAD);
			count++;
			seqs.push_back(sub);
			if(count % 10000 == 0)
			{
				printf("sent %d sequences\n", count);
			}
		}
	}
	
	commLayer.SendControlMessage(1, APM_DONELOAD);
	printf("seq load done, %d sequences parsed\n", count);
	
	count = 0;
	for(PSequenceVector::iterator iter = seqs.begin(); iter != seqs.end(); iter++)
	{
		bool exists = commLayer.CheckForSequence(1, *iter);
		assert(exists);
		
		if(count % 10000 == 0)
		{
			printf("passed %d sequences\n", count);
		}
		count++;
	}
	
	commLayer.SendControlMessage(1, APM_FINISHED);
}

void printUsage()
{
	printf("usage: ABYSS <reads fasta file> <kmer size> <number of trimming steps>\n");	
}
