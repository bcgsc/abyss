#include <stdio.h>

#include <vector>
#include <stdio.h>
#include <mpi.h>
#include "parallelAbyss.h"
#include "CommonUtils.h"
#include "FastaReader.h"
#include "NetworkSequenceCollection.h"

int rank;

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
	//int numTrims = atoi(argv[3]);
	
	
   double starttime, endtime; 
   starttime = MPI::Wtime();
	
	// start mpi process
	MPI::Init(argc,argv);
	
	// get my rank
	rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	NetworkSequenceCollection networkSeqs(rank, size, kmerSize, readLen);
	if(rank == 0)
	{
		printf("num nodes: %d\n", size);
	}
	
	if(rank == 0)
	{
		networkSeqs.runControl(fastaFile, readLen, kmerSize);
	}
	else
	{
		printf("%d running\n", rank);
		networkSeqs.run(readLen, kmerSize);
	}
	printf("%d exitting\n", rank);
	
	endtime = MPI::Wtime();
	
	printf("%d elapsed seconds: %f\n", rank, endtime - starttime);
	MPI::Finalize();
	return 0;
}

void printUsage()
{
	printf("usage: ABYSS <reads fasta file> <kmer size> <number of trimming steps>\n");	
}
