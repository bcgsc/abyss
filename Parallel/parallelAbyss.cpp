#include <cstdio>
#include <sstream>
#include <vector>
#include <mpi.h>
#include "parallelAbyss.h"
#include "CommonUtils.h"
#include "FastaReader.h"
#include "NetworkSequenceCollection.h"
#include "Options.h"
#include "ParallelFastaWriter.h"
#include "Timer.h"

/** MPI size */
static int mpi_size;

static void concatenateContigs()
{
	puts("Concatenating contig files");
	int ret = system("echo -n >pcontigs.fa");
	assert(ret == 0);
	(void)ret;
	for (int i = 0; i < mpi_size; i++) {
		std::ostringstream s;
		s << "cat contigs-" << i << ".fa >>pcontigs.fa";
		ret = system(s.str().c_str());
		assert(ret == 0);
	}
}

int main(int argc, char** argv)
{	
	Timer timer("ParallelAbyss");
	
	// Set stdout to be line buffered.
	setvbuf(stdout, NULL, _IOLBF, 0);

	// start mpi process
	MPI_Init(&argc,&argv);
	
	// get my rank and the world size
	MPI_Comm_rank(MPI_COMM_WORLD, &opt::rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	
	opt::parse(argc, argv);
	if (opt::rank == 0)
		printf("Running on %d processors\n", mpi_size);

	NetworkSequenceCollection networkSeqs(opt::rank, mpi_size,
			opt::kmerSize, opt::readLen);

	if (opt::rank == 0) {
		networkSeqs.runControl();
		concatenateContigs();
		puts("Done.");
	} else
		networkSeqs.run();

	MPI_Finalize();
	return 0;
}
