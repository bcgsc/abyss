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
#include "Timer.h"

int rank;

int main(int argc, char** argv)
{	
	opt::parse(argc, argv);

	Timer timer("ParallelAbyss");
	
	// start mpi process
	MPI_Init(&argc,&argv);
	
	// get my rank and the world size
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	

	NetworkSequenceCollection networkSeqs(rank, size,
			opt::kmerSize, opt::readLen);

	if (rank == 0) {
		printf("num nodes: %d\n", size);
		networkSeqs.runControl();
	} else
		networkSeqs.run();

	MPI_Finalize();
	return 0;
}
