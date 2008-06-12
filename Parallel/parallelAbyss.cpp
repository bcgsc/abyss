#include <cstdio>
#include <vector>
#include <mpi.h>
#include "parallelAbyss.h"
#include "CommonUtils.h"
#include "FastaReader.h"
#include "NetworkSequenceCollection.h"
#include "Options.h"
#include "ParallelFastaWriter.h"
#include "Timer.h"

int main(int argc, char** argv)
{	
	Timer timer("ParallelAbyss");
	
	// Set stdout to be line buffered.
	setvbuf(stdout, NULL, _IOLBF, 0);

	// start mpi process
	MPI_Init(&argc,&argv);
	
	// get my rank and the world size
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &opt::rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	opt::parse(argc, argv);
	if (opt::rank == 0)
		printf("Running on %d processors\n", size);
	if (opt::snpFile != NULL)
		freopen(NULL, "a", opt::snpFile);

	NetworkSequenceCollection networkSeqs(opt::rank, size,
			opt::kmerSize, opt::readLen);

	if (opt::rank == 0)
		networkSeqs.runControl();
	else
		networkSeqs.run();

	MPI_Finalize();
	return 0;
}
