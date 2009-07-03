#include "Aligner.h"
#include "PairedAlgorithms.h"
#include "PrefixIterator.h"
#include "FastaReader.h"
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <string>
#include <pthread.h>
#include <semaphore.h>

using namespace std;

#define PROGRAM "KAligner"

static const char *VERSION_MESSAGE =
PROGRAM " (ABySS) " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2009 Canada's Michael Smith Genome Science Centre\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION]... QUERY TARGET\n"
"Align the sequences of QUERY against those of TARGET.\n"
"All perfect matches of at least k bases will be found.\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -j, --threads=THREADS the max number of threads created\n"
"                        set to 0 for one thread per reads file\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static unsigned k;
	static int threads = 1;
	static int verbose;
	extern bool colourSpace;
}

static const char* shortopts = "k:o:j:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "threads",     required_argument,	NULL, 'j' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

// Functions
static void readContigsIntoDB(string refFastaFile, Aligner& aligner);
void *alignReadsToDB(void *arg);

// Global variables
static Aligner *aligner;
static unsigned readCount;
static pthread_mutex_t mutexCout, mutexCerr;
static pthread_attr_t attr;
static sem_t activeThreads;

static void getReadFiles(string readsFile, vector<pthread_t>* threads)
{
	// Ensure we don't create more than opt::threads threads at a time.
	if (opt::threads > 0)
		sem_wait(&activeThreads);

	if (opt::verbose > 0) {
		pthread_mutex_lock(&mutexCerr);
		cerr << "Reading `" << readsFile << "'...\n";
		pthread_mutex_unlock(&mutexCerr);
	}

	// Create a copy of the readsFile string to prevent it from being
	// overwriten in next iteration.
	string* temp = new string(readsFile);

	pthread_t thread;
	pthread_create(&thread, &attr, alignReadsToDB, (void*) temp);
	threads->push_back(thread);
}

int main(int argc, char** argv)
{
	bool die = false;
	for (char c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
			case 'j': arg >> opt::threads; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (opt::k <= 0) {
		cerr << PROGRAM ": missing -k,--kmer option\n";
		die = true;
	}

	if (argc - optind < 2) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	string refFastaFile(argv[argc - 1]);

	if (opt::verbose > 0)
		cerr << "k: " << opt::k
			<< " Target: " << refFastaFile
			<< endl;

	aligner = new Aligner(opt::k);
	readContigsIntoDB(refFastaFile, *aligner);

	// Need to initialize mutex's before threads are created.
	pthread_mutex_init(&mutexCout, NULL);
	pthread_mutex_init(&mutexCerr, NULL);
	if (opt::threads > 0)
		sem_init(&activeThreads, 0, opt::threads);

	// Make the threads joinable so this main thread can wait untill
	// they are finished.
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	readCount = 0;
	vector<pthread_t> threads;
	for_each(argv + optind, argv + argc - 1,
			bind2nd(ptr_fun(getReadFiles), &threads));

	void *status;
	// Wait for all threads to finish.
	for (size_t i = 0; i < threads.size(); i++)
		pthread_join(threads[i], &status);

	if (opt::verbose > 0)
		cerr << "Aligned " << readCount << " reads\n";

	pthread_exit(NULL);
}

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

static void printProgress(const Aligner& align, unsigned count)
{
	size_t size = align.size();
	size_t buckets = align.bucket_count();
	cerr << "Read " << count << " contigs. "
		<< "Hash load: " << size <<
		" / " << buckets << " = " << (float)size / buckets << endl;
}

static void readContigsIntoDB(string refFastaFile, Aligner& aligner)
{
	int count = 0;
	ifstream fileHandle(refFastaFile.c_str());
	assert_open(fileHandle, refFastaFile);

	while(!fileHandle.eof() && fileHandle.peek() != EOF)
	{
		ContigID contigID;
		Sequence seq;
		int length;
		double coverage;

		PairedAlgorithms::parseContigFromFile(fileHandle, contigID, seq, length, coverage);

		if (count == 0) {
			// Detect colour-space contigs.
			opt::colourSpace = isdigit(seq[0]);
		} else {
			if (opt::colourSpace)
				assert(isdigit(seq[0]));
			else
				assert(isalpha(seq[0]));
		}

		aligner.addReferenceSequence(contigID, seq);

		count++;
		if (opt::verbose > 0 && count % 100000 == 0)
			printProgress(aligner, count);
	}
	if (opt::verbose > 0)
			printProgress(aligner, count);

	fileHandle.close();
}

void *alignReadsToDB(void *readsFile)
{
	FastaReader fileHandle(((string*)readsFile)->c_str());

	ostringstream output;
	prefix_ostream_iterator<Alignment> out(output, "\t");

	while (fileHandle.isGood()) {
		string id;
		Sequence seq = fileHandle.ReadSequence(id);

		if (opt::colourSpace)
			assert(isdigit(seq[0]));
		else
			assert(isalpha(seq[0]));

		size_t pos = seq.find_first_not_of("ACGT0123");
		if (pos == string::npos)
			aligner->alignRead(seq, out);

		pthread_mutex_lock(&mutexCout);
		cout << id << output.str() << '\n';
		assert(cout.good());
		pthread_mutex_unlock(&mutexCout);

		output.str("");

		if (opt::verbose > 0) {
			pthread_mutex_lock(&mutexCerr);
			if (++readCount % 1000000 == 0)
				cerr << "Aligned " << readCount << " reads\n";
			pthread_mutex_unlock(&mutexCerr);
		}
	}
	if (opt::threads > 0)
		sem_post(&activeThreads);
	pthread_exit(NULL);
}
