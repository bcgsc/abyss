#include "Aligner.h"
#include "Align/Options.h"
#include "Common/Options.h"
#include "DataLayer/Options.h"
#include "AffixIterator.h"
#include "Barrier.h"
#include "FastaReader.h"
#include "Uncompress.h"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <pthread.h>
#include <semaphore.h>
#include <sstream>
#include <string>
#include <sys/stat.h>

using namespace std;

#define PROGRAM "KAligner"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2010 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... QUERY... TARGET\n"
"Align the sequences of the files QUERY against those of the file TARGET.\n"
"All perfect matches of at least k bases will be found.\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"      --no-multimap     disallow duplicate k-mer in the target [default]\n"
"  -m, --multimap        allow duplicate k-mer in the target\n"
"  -i, --ignore-multimap ignore duplicate k-mer in the target\n"
"  -j, --threads=N       use N threads [2] up to one per query file\n"
"                        or if N is 0 use one thread per query file\n"
#if _POSIX_BARRIERS > 0
"      --sync=COUNT      synchronize threads every COUNT alignments [10000]\n"
"      --no-sync         do not synchronize threads\n"
#endif
"  -v, --verbose         display verbose output\n"
"      --seq             print the sequence with the alignments\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static unsigned k;
	static int threads = 2;
	static int printSeq;

	/** Synchronize the threads with a barrier. */
	static int sync = 10000;
}

static const char shortopts[] = "ik:mo:j:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_SYNC };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "no-multi",    no_argument,     &opt::multimap, opt::ERROR },
	{ "multimap",    no_argument,     &opt::multimap, opt::MULTIMAP },
	{ "ignore-multimap", no_argument, &opt::multimap, opt::IGNORE },
	{ "sync",        required_argument, NULL, OPT_SYNC },
	{ "no-sync",     no_argument,       &opt::sync, 0 },
	{ "threads",     required_argument,	NULL, 'j' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "seq",		 no_argument,		&opt::printSeq, 1 },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Return the number of k-mer in the specified file. */
static size_t countKmer(const string& path)
{
	struct stat st;
	if (stat(path.c_str(), &st) == -1) {
		perror(path.c_str());
		exit(EXIT_FAILURE);
	}

	if (!S_ISREG(st.st_mode)) {
		cerr << "Not calculating k-mer in `" << path
			<< "', because it is not a regular file.\n";
		return 500000000;
	}

	if (opt::verbose > 0)
		cerr << "Reading target `" << path << "'..." << endl;

	ifstream in(path.c_str());
	assert(in.is_open());
	size_t scaffolds = 0, contigs = 0, bases = 0;
	enum { ID, SEQUENCE, GAP } state = SEQUENCE;
	for (char c; in.get(c);) {
		c = toupper(c);
		switch (state) {
		  case ID:
			if (c == '\n')
				state = SEQUENCE;
			break;
		  case SEQUENCE:
		  case GAP:
			switch (c) {
			  case '>':
				scaffolds++;
				contigs++;
				state = ID;
				break;
			  case 'N':
				if (state != GAP)
					contigs++;
				state = GAP;
				break;
			  case 'A': case 'C': case 'G': case 'T':
				bases++;
				state = SEQUENCE;
				break;
			}
			break;
		}
	}

	size_t overlaps = contigs * (opt::k-1);
	size_t kmer = bases - overlaps;
	if (opt::verbose > 0)
		cerr << "Read " << bases << " bases, "
			<< contigs << " contigs, " << scaffolds << " scaffolds"
			" from `" << path << "'. "
			"Expecting " << kmer << " k-mer.\n";
	assert(bases > overlaps);
	return kmer;
}

template <class SeqPosHashMap>
static void readContigsIntoDB(string refFastaFile,
		Aligner<SeqPosHashMap>& aligner);
void *alignReadsToDB(void *arg);

/** Unique aligner using map */
static Aligner<SeqPosHashUniqueMap> *g_aligner_u;

/** Multimap aligner using multimap */
static Aligner<SeqPosHashMultiMap> *g_aligner_m;

/** Number of reads. */
static unsigned g_readCount;

/** Number of reads that aligned. */
static unsigned g_alignedCount;

static pthread_mutex_t g_mutexCout, g_mutexCerr;
static sem_t g_activeThreads;

static pthread_t getReadFiles(const char *readsFile)
{
	// Ensure we don't create more than opt::threads threads at a time.
	sem_wait(&g_activeThreads);

	if (opt::verbose > 0) {
		pthread_mutex_lock(&g_mutexCerr);
		cerr << "Reading `" << readsFile << "'...\n";
		pthread_mutex_unlock(&g_mutexCerr);
	}

	pthread_t thread;
	pthread_create(&thread, NULL, alignReadsToDB, (void*)readsFile);
	return thread;
}

/** A barrier used to synchronize the worker threads. */
static Barrier g_barrier;

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
			case 'm': opt::multimap = opt::MULTIMAP; break;
			case 'i': opt::multimap = opt::IGNORE; break;
			case 'j': arg >> opt::threads; break;
			case 'v': opt::verbose++; break;
			case OPT_SYNC: arg >> opt::sync; break;
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
	Kmer::setLength(opt::k);

	if (argc - optind < 2) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	string refFastaFile(argv[--argc]);

	int numQuery = argc - optind;
	if (opt::threads <= 0 || opt::threads > numQuery)
		opt::threads = numQuery;
	if (opt::threads == 1)
		opt::sync = 0;
#if !_POSIX_BARRIERS
	opt::sync = 0;
#endif

	size_t numKmer = countKmer(refFastaFile);
	if (opt::multimap == opt::MULTIMAP) {
		g_aligner_m = new Aligner<SeqPosHashMultiMap>(opt::k,
				numKmer);
		readContigsIntoDB(refFastaFile, *g_aligner_m);
	} else {
#if HAVE_GOOGLE_SPARSE_HASH_SET
		g_aligner_u = new Aligner<SeqPosHashUniqueMap>(opt::k,
				numKmer, 0.2);
#else
		g_aligner_u = new Aligner<SeqPosHashUniqueMap>(opt::k,
				numKmer);
#endif
		readContigsIntoDB(refFastaFile, *g_aligner_u);
	}

	// Need to initialize mutex's before threads are created.
	pthread_mutex_init(&g_mutexCout, NULL);
	pthread_mutex_init(&g_mutexCerr, NULL);
	sem_init(&g_activeThreads, 0, opt::threads);
	if (opt::sync) {
		// Initialize the barrier for the initial number of threads,
		// then set the thread count to 0. Each thread increments the
		// barrier upon its creation and decrements the barrier upon
		// its completion.
		g_barrier.init(opt::threads);
		g_barrier = 0;
	}

	g_readCount = 0;
	vector<pthread_t> threads;
	transform(argv + optind, argv + argc, back_inserter(threads),
			getReadFiles);

	if (opt::sync) {
		for (int i = 0; i < opt::threads; i++) {
			sem_wait(&g_activeThreads);
			// Signal the barrier to indicate that the worker thread
			// has completed and there won't be another to replace it.
			g_barrier.signal();
		}
	}

	void *status;
	// Wait for all threads to finish.
	for (size_t i = 0; i < threads.size(); i++)
		pthread_join(threads[i], &status);

	if (opt::verbose > 0)
		cerr << "Aligned " << g_alignedCount
			<< " of " << g_readCount << " reads ("
			<< (float)100 * g_alignedCount / g_readCount << "%)\n";

	if (opt::multimap == opt::MULTIMAP)
		delete g_aligner_m;
	else
		delete g_aligner_u;

	return 0;
}

template <class SeqPosHashMap>
static void printProgress(const Aligner<SeqPosHashMap>& align,
		unsigned count)
{
	size_t size = align.size();
	size_t buckets = align.bucket_count();
	cerr << "Read " << count << " contigs. "
		<< "Hash load: " << size <<
		" / " << buckets << " = " << (float)size / buckets << endl;
}

template <class SeqPosHashMap>
static void readContigsIntoDB(string refFastaFile,
		Aligner<SeqPosHashMap>& aligner)
{
	if (opt::verbose > 0)
		cerr << "Reading target `" << refFastaFile << "'..." << endl;

	unsigned count = 0;
	FastaReader in(refFastaFile.c_str(),
			FastaReader::KEEP_N | FastaReader::FOLD_CASE);
	for (FastaRecord rec; in >> rec;) {
		if (count == 0) {
			// Detect colour-space contigs.
			opt::colourSpace = isdigit(rec.seq[0]);
		} else {
			if (opt::colourSpace)
				assert(isdigit(rec.seq[0]));
			else
				assert(isalpha(rec.seq[0]));
		}

		aligner.addReferenceSequence(rec.id, rec.seq);

		count++;
		if (opt::verbose > 0 && count % 100000 == 0)
			printProgress(aligner, count);
	}
	assert(in.eof());
	if (opt::verbose > 0)
		printProgress(aligner, count);
}

void *alignReadsToDB(void* readsFile)
{
	int barrierCount = 0;
	if (opt::sync)
		++g_barrier;

	opt::chastityFilter = false;
	opt::trimMasked = false;

	// Lock `uncompress', which is not thread safe.
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_lock(&mutex);
	FastaReader fileHandle((const char *)readsFile,
			FastaReader::KEEP_N | FastaReader::FOLD_CASE);
	pthread_mutex_unlock(&mutex);

	for (FastaRecord rec; fileHandle >> rec;) {
		const Sequence& seq = rec.seq;
		ostringstream output;
		if (seq.find_first_not_of("ACGT0123") == string::npos) {
			if (opt::colourSpace)
				assert(isdigit(seq[0]));
			else
				assert(isalpha(seq[0]));

			if (opt::multimap == opt::MULTIMAP)
				g_aligner_m->alignRead(seq,
						affix_ostream_iterator<Alignment>(
							output, "\t"));
			else
				g_aligner_u->alignRead(seq,
						affix_ostream_iterator<Alignment>(
							output, "\t"));
		}

		pthread_mutex_lock(&g_mutexCout);
		cout << rec.id;
		if (opt::printSeq) {
			cout << ' ';
			if (opt::colourSpace)
				cout << rec.anchor;
			cout << seq;
		}
		string s = output.str();
		cout << s << '\n';
		assert(cout.good());
		pthread_mutex_unlock(&g_mutexCout);

		if (opt::verbose > 0) {
			pthread_mutex_lock(&g_mutexCerr);
			if (!s.empty())
				g_alignedCount++;
			if (++g_readCount % 1000000 == 0)
				cerr << "Aligned " << g_readCount << " reads\n";
			pthread_mutex_unlock(&g_mutexCerr);
		}

		// Synchronize the threads periodically.
		if (opt::sync && ++barrierCount == opt::sync) {
			barrierCount = 0;
			g_barrier.wait();
		}
	}
	assert(fileHandle.eof());
	if (opt::sync)
		--g_barrier;
	sem_post(&g_activeThreads);
	return NULL;
}
