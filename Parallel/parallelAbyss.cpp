#include "config.h"
#include <climits> // for HOST_NAME_MAX
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <unistd.h> // for gethostname
#include <vector>
#include <mpi.h>
#include "FastaReader.h"
#include "Log.h"
#include "NetworkSequenceCollection.h"
#include "Options.h"
#include "Timer.h"

using namespace std;

static void concatenateFiles(string dest,
		string prefix, string suffix)
{
	printf("Concatenating to %s\n", dest.c_str());
	std::ostringstream s;
	s << "cat";
	for (int i = 0; i < opt::numProc; i++)
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
	Timer timer("Total");

	// Set stdout to be line buffered.
	setvbuf(stdout, NULL, _IOLBF, 0);

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &opt::rank);
	MPI_Comm_size(MPI_COMM_WORLD, &opt::numProc);

	opt::parse(argc, argv);
	if (opt::rank == 0)
		printf("Running on %d processors\n", opt::numProc);

	MPI_Barrier(MPI_COMM_WORLD);
	char hostname[HOST_NAME_MAX];
	gethostname(hostname, sizeof hostname);
	PrintDebug(0, "Running on host %s\n", hostname);
	MPI_Barrier(MPI_COMM_WORLD);

	NetworkSequenceCollection networkSeqs(opt::rank, opt::numProc);

	if (opt::rank == 0)
		networkSeqs.runControl();
	else
		networkSeqs.run();

	if (opt::rank == 0) {
		concatenateFiles(opt::contigsPath, "contigs-", ".fa");
		if (opt::snpPath.length() > 0)
			concatenateFiles(opt::snpPath, "snp-", ".fa");
		puts("Done.");
	}
	
	MPI_Finalize();

	return 0;
}
