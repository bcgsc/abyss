#include "config.h"
#include "Log.h"
#include "NetworkSequenceCollection.h"
#include "Assembly/Options.h"
#include "Common/Options.h"
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
#include <unistd.h> // for gethostname and sync
#include <vector>

using namespace std;

/** Execute a shell command and check its return status. */
static void systemx(const string& command)
{
	if (opt::verbose > 0)
		cout << command << endl;
	int ret = system(command.c_str());
	if (ret == 0)
		return;
	cerr << "error: command failed: `" << command << "'\n";
	if (ret == -1)
		cerr << "system() failed: " << strerror(errno) << endl;
	exit(ret == -1 ? EXIT_FAILURE : ret);
}

/** Concatenate files using the specified command and remove them. */
static void concatenateFiles(const string& dest,
		const string& prefix, const string& suffix,
		const string& command = "cat")
{
	cout << "Concatenating to " << dest << endl;
	ostringstream s;
	s << command;
	for (int i = 0; i < opt::numProc; i++)
		s << ' ' << prefix << i << suffix;
	s << " >'" << dest << '\'';
	systemx(s.str());

	s.str("");
	s << "rm";
	for (int i = 0; i < opt::numProc; i++)
		s << ' ' << prefix << i << suffix;
	systemx(s.str());
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
		cout << "Running on " << opt::numProc << " processors\n";

	MPI_Barrier(MPI_COMM_WORLD);
	char hostname[HOST_NAME_MAX];
	gethostname(hostname, sizeof hostname);
	logger(0) << "Running on host " << hostname << endl;
	MPI_Barrier(MPI_COMM_WORLD);

	NetworkSequenceCollection networkSeqs;
	if (opt::rank == 0)
		networkSeqs.runControl();
	else
		networkSeqs.run();

	cout << "Synchronizing file system...\n;
	sync();
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();

	if (opt::rank == 0) {
		concatenateFiles(opt::contigsPath, "contigs-", ".fa",
				"awk '/^>/ { $1=\">\" i++ } { print }'");
		if (!opt::snpPath.empty())
			concatenateFiles(opt::snpPath, "snp-", ".fa");
		cout << "Done." << endl;
	}

	return 0;
}
