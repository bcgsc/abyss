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
	Timer timer("ParallelAbyss");
	
	// start mpi process
	MPI_Init(&argc,&argv);
	
	// get my rank and the world size
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	if (rank == 0) {
		for (int i = 0; i < argc; i++) {
			if (i != 0)
				putchar(' ');
			fputs(argv[i], stdout);
		}
		putchar('\n');
		printf("Running on %d processors\n", size);
	}

	opt::parse(argc, argv);
	NetworkSequenceCollection networkSeqs(rank, size,
			opt::kmerSize, opt::readLen);

	if (rank == 0)
		networkSeqs.runControl();
	else
		networkSeqs.run();

	MPI_Finalize();
	return 0;
}
