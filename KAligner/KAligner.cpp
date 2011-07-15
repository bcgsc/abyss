#include "Aligner.h"
#include "Align/Options.h"
#include "Barrier.h"
#include "Common/Options.h"
#include "DataLayer/Options.h"
#include "FastaReader.h"
#include "Iterator.h"
#include "SAM.h"
#include "StringUtil.h" // for toSI
#include "Uncompress.h"
#include "Pipe.h"
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
#include <unistd.h> // for sbrk

using namespace std;

#define PROGRAM "KAligner"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2011 Canada's Michael Smith Genome Science Centre\n";

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
"      --no-sam          output the results in KAligner format\n"
"      --sam             output the results in SAM format\n"
"      --seq             print the sequence with the alignments\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

/** Enumeration of output formats */
enum format { KALIGNER, SAM };

namespace opt {
	static unsigned k;
	static int threads = 2;
	static int printSeq;

	/** Synchronize the threads with a barrier. */
	static int sync = 10000;

	/** Output formats */
	static int format;
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
	{ "no-sam",      no_argument,       &opt::format, KALIGNER },
	{ "sam",         no_argument,       &opt::format, SAM },
	{ "no-seq",		 no_argument,		&opt::printSeq, 0 },
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
			  case 'B': case 'D': case 'H': case 'K': case 'M':
			  case 'R': case 'S': case 'V': case 'W': case 'Y':
				if (state != GAP)
					contigs++;
				state = GAP;
				break;
			  case 'A': case 'C': case 'G': case 'T':
			  case '0': case '1': case '2': case '3':
				bases++;
				state = SEQUENCE;
				break;
			  case '\n':
				break;
			  default:
				cerr << "error: unexpected character: `" << c << "'\n";
				exit(EXIT_FAILURE);
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
static void *alignReadsToDB(void *arg);

/** Unique aligner using map */
static Aligner<SeqPosHashUniqueMap> *g_aligner_u;

/** Multimap aligner using multimap */
static Aligner<SeqPosHashMultiMap> *g_aligner_m;

/** Number of reads. */
static unsigned g_readCount;

/** Number of reads that aligned. */
static unsigned g_alignedCount;

static pthread_mutex_t g_mutexIn, g_mutexCout, g_mutexCerr;
static vector<const char *> g_files;
static FastaReader* g_fileHandle;

static void getReadFiles(const char *readsFile)
{
	if (opt::verbose > 0) {
		pthread_mutex_lock(&g_mutexCerr);
		cerr << "Reading `" << readsFile << "'...\n";
		pthread_mutex_unlock(&g_mutexCerr);
	}

	if (g_fileHandle != NULL) {
		delete g_fileHandle;
		g_fileHandle = NULL;
	}

	assert(g_fileHandle == NULL);

	g_fileHandle = new FastaReader(readsFile,
			FastaReader::FOLD_CASE);

	//pthread_t thread;
	//pthread_create(&thread, NULL, alignReadsToDB, (void*)readsFile);
	//return p_fr;
}

int main(int argc, char** argv)
{
	string commandLine;
	{
		ostringstream ss;
		char** last = argv + argc - 1;
		copy(argv, last, ostream_iterator<const char *>(ss, " "));
		ss << *last;
		commandLine = ss.str();
	}

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
	if (opt::threads <= 0)
		opt::threads = numQuery;
	if (opt::threads == 1)
		opt::sync = 0;
#if !_POSIX_BARRIERS
	opt::sync = 0;
#endif

	// SAM headers.
	cout << "@HD\tVN:1.0\n"
		"@PG\tID:" PROGRAM "\tVN:" VERSION "\t"
		"CL:" << commandLine << '\n';

	size_t numKmer = countKmer(refFastaFile);
	if (opt::multimap == opt::MULTIMAP) {
		g_aligner_m = new Aligner<SeqPosHashMultiMap>(opt::k,
				numKmer);
		readContigsIntoDB(refFastaFile, *g_aligner_m);
	} else {
#if HAVE_GOOGLE_SPARSE_HASH_MAP
		g_aligner_u = new Aligner<SeqPosHashUniqueMap>(opt::k,
				numKmer, 0.3);
#else
		g_aligner_u = new Aligner<SeqPosHashUniqueMap>(opt::k,
				numKmer);
#endif
		readContigsIntoDB(refFastaFile, *g_aligner_u);
	}

	// Need to initialize mutex's before threads are created.
	pthread_mutex_init(&g_mutexCout, NULL);
	pthread_mutex_init(&g_mutexCerr, NULL);
	pthread_mutex_init(&g_mutexIn, NULL);

	g_readCount = 0;
	for (int i=optind; i < argc; i++)
		g_files.push_back(argv[i]);

	getReadFiles(g_files[0]);
	g_files.erase(g_files.begin());
	
	vector<pthread_t> threads;
	for (int i = 0; i < opt::threads; i++) {
		pthread_t thread;
		pthread_create(&thread, NULL, alignReadsToDB, NULL);
		threads.push_back(thread);
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
	delete g_fileHandle;

	assert(g_files.empty());
	return 0;
}

/** Start of the data segment. */
static intptr_t sbrk0 = reinterpret_cast<intptr_t>(sbrk(0));

template <class SeqPosHashMap>
static void printProgress(const Aligner<SeqPosHashMap>& align,
		unsigned count)
{
	ptrdiff_t bytes = reinterpret_cast<intptr_t>(sbrk(0)) - sbrk0;
	size_t size = align.size();
	size_t buckets = align.bucket_count();
	cerr << "Read " << count << " contigs. "
		"Hash load: " << size << " / " << buckets
		<< " = " << (float)size / buckets
		<< " using " << toSI(bytes) << "B." << endl;
}

template <class SeqPosHashMap>
static void readContigsIntoDB(string refFastaFile,
		Aligner<SeqPosHashMap>& aligner)
{
	if (opt::verbose > 0)
		cerr << "Reading target `" << refFastaFile << "'..." << endl;

	unsigned count = 0;
	FastaReader in(refFastaFile.c_str(), FastaReader::FOLD_CASE);
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

		cout << "@SQ\tSN:" << rec.id << "\tLN:" << rec.seq.length() << '\n';
		aligner.addReferenceSequence(rec.id, rec.seq);

		count++;
		if (opt::verbose > 0 && count % 100000 == 0)
			printProgress(aligner, count);
	}
	assert(in.eof());
	if (opt::verbose > 0)
		printProgress(aligner, count);

	if (opt::multimap == opt::IGNORE) {
		// Count the number of duplicate k-mer in the target.
		size_t duplicates = aligner.countDuplicates();
		if (duplicates > 0)
			cerr << "Found " << duplicates
				<< " (" << (float)100 * duplicates / aligner.size()
				<< "%) duplicate k-mer.\n";
	}
}

static bool getNextRec(FastaRecord& rec)
{
	pthread_mutex_lock(&g_mutexIn);
	FastaReader& fileHandle = *g_fileHandle;
	while (!(fileHandle >> rec)) {
		assert(fileHandle.eof());
		if (g_files.empty()) {
			pthread_mutex_unlock(&g_mutexIn);
			return false;
		}
		getReadFiles(g_files[0]);
		g_files.erase(g_files.begin());
	}

	pthread_mutex_unlock(&g_mutexIn);
	return true;
}

static void *alignReadsToDB(void*)
{
	opt::chastityFilter = false;
	opt::trimMasked = false;

	for (FastaRecord rec; getNextRec(rec);) {
		const Sequence& seq = rec.seq;
		ostringstream output;
		if (seq.find_first_not_of("ACGT0123") == string::npos) {
			if (opt::colourSpace)
				assert(isdigit(seq[0]));
			else
				assert(isalpha(seq[0]));
		}

		switch (opt::format) {
		  case KALIGNER:
			if (opt::multimap == opt::MULTIMAP)
				g_aligner_m->alignRead(rec.id, seq,
						affix_ostream_iterator<Alignment>(
							output, "\t"));
			else
				g_aligner_u->alignRead(rec.id, seq,
						affix_ostream_iterator<Alignment>(
							output, "\t"));
			break;
		  case SAM:
			if (opt::multimap == opt::MULTIMAP)
				g_aligner_m->alignRead(rec.id, seq,
						ostream_iterator<SAMRecord>(output, "\n"));
			else
				g_aligner_u->alignRead(rec.id, seq,
						ostream_iterator<SAMRecord>(output, "\n"));
			break;
		}

		string s = output.str();
		pthread_mutex_lock(&g_mutexCout);
		switch (opt::format) {
		  case KALIGNER:
			cout << rec.id;
			if (opt::printSeq) {
				cout << ' ';
				if (opt::colourSpace)
					cout << rec.anchor;
				cout << seq;
			}
			cout << s << '\n';
			break;
		  case SAM:
			cout << s;
			break;
		}
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

	}
	return NULL;
}
