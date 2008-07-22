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

static void concatenateFiles(string dest,
		string prefix, string suffix)
{
	printf("Concatenating to %s\n", dest.c_str());
	std::ostringstream s;
	s << "cat";
	for (int i = 0; i < mpi_size; i++)
		s << ' ' << prefix << i << suffix;
	s << " >'" << dest << '\'';
	if (opt::verbose > 0)
		puts(s.str().c_str());
	int ret = system(s.str().c_str());
	if (ret != 0) {
		fprintf(stderr, "error: command failed: %s\n",
				s.str().c_str());
		if (ret == -1)
			perror("system");
		exit(ret == -1 ? EXIT_FAILURE : ret);
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

	if (opt::rank == 0)
		networkSeqs.runControl();
	else
		networkSeqs.run();

	MPI_Finalize();

	if (opt::rank == 0) {
		concatenateFiles("pcontigs.fa", "contigs-", ".fa");
		if (opt::snpPath.length() > 0)
			concatenateFiles(opt::snpPath, "snp-", ".fa");
		puts("Done.");
	}

	return 0;
}
