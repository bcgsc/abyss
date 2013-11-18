#include "config.h"
#include "Log.h"
#include "NetworkSequenceCollection.h"
#include "Assembly/Options.h"
#include "Common/Options.h"
#include "DataLayer/FastaReader.h"
#include "Timer.h"
#include "Uncompress.h"
#include <cerrno>
#include <climits> // for HOST_NAME_MAX
#include <cstdio> // for setvbuf
#include <cstdlib>
#include <cstring> // for strerror
#include <iostream>
#include <mpi.h>
#include <sstream>
#include <unistd.h> // for gethostname
#include <vector>

using namespace std;

static const char* FASTA_SUFFIX = ".fa";

static void mergeFastaFiles(const string& outputPath, const string& inputPathPrefix, bool generateNewIds = false)
{
	cout << "Concatenating fasta files to " << outputPath << endl;

	// write merged FASTA file

	FastaWriter writer(outputPath.c_str());
	uint64_t seqid = 0;
	for(int i = 0; i < opt::numProc; i++) {
		ostringstream filename;
		filename << inputPathPrefix << i << FASTA_SUFFIX;
		assert(filename.good());
		FastaReader reader(filename.str().c_str(), FastaReader::NO_FOLD_CASE);
		for (FastaRecord rec; reader >> rec;) {
			if (generateNewIds)
				writer.WriteSequence(rec.seq, seqid++, rec.comment);
			else
				writer.WriteSequence(rec.seq, rec.id, rec.comment);
		}
		assert(reader.eof());
	}

	// remove temp FASTA files

	bool die = false;
	for (int i = 0; i < opt::numProc; i++) {
		ostringstream s;
		s << inputPathPrefix << i << FASTA_SUFFIX;
		const char* path = s.str().c_str();
		if (unlink(path) == -1) {
			cerr << "error: removing `" << path << "': "
				<< strerror(errno) << endl;
			die = true;
		}
	}
	if (die)
		exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
	Timer timer("Total");

	// Set stdout to be line buffered.
	setvbuf(stdout, NULL, _IOLBF, 0);

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &opt::rank);
	MPI_Comm_size(MPI_COMM_WORLD, &opt::numProc);

	// OMPI-1.6.1 and later reset the SIGCHLD handler so we need to
	// reinitialize uncompress.
	uncompress_init();

	opt::parse(argc, argv);
	if (opt::rank == 0)
		cout << "Running on " << opt::numProc << " processors\n";

	MPI_Barrier(MPI_COMM_WORLD);
	char hostname[HOST_NAME_MAX];
	gethostname(hostname, sizeof hostname);
	logger(0) << "Running on host " << hostname << endl;
	MPI_Barrier(MPI_COMM_WORLD);

	if (opt::rank == 0) {
		NetworkSequenceCollection networkSeqs;
		networkSeqs.runControl();
	} else {
		NetworkSequenceCollection networkSeqs;
		networkSeqs.run();
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	if (opt::rank == 0) {
		mergeFastaFiles(opt::contigsPath, "contigs-", true);
		if (!opt::snpPath.empty())
			mergeFastaFiles(opt::snpPath, "snp-");
		cout << "Done." << endl;
	}

	return 0;
}
