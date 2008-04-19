#include <stdio.h>

#include <vector>
#include <stdio.h>
#include <mpi.h>
#include "parallelAbyss.h"
#include "CommonUtils.h"
#include "FastaReader.h"
#include "NetworkSequenceCollection.h"
#include "Options.h"
#include "ParallelFastaWriter.h"

int rank;

int main(int argc, char** argv)
{	
	opt::parse(argc, argv);

   double starttime, endtime; 
   starttime = MPI::Wtime();
	
	// start mpi process
	MPI::Init(argc,argv);
	
	// get my rank
	rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	NetworkSequenceCollection networkSeqs(rank, size,
			opt::kmerSize, opt::readLen);
	if(rank == 0)
	{
		printf("num nodes: %d\n", size);
	}
	
	if(rank == 0)
	{
		networkSeqs.runControl(opt::fastaFile,
				opt::readLen, opt::kmerSize);
	}
	else
	{
		networkSeqs.run(opt::readLen, opt::kmerSize);
	}

	endtime = MPI::Wtime();
	
	printf("%d elapsed seconds: %f\n", rank, endtime - starttime);
	MPI::Finalize();
	return 0;
}
