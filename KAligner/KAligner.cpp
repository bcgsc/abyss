#include "Aligner.h"
#include "Common/Options.h"
#include "DataLayer/Options.h"
#include "KAligner/Options.h"
#include "FastaReader.h"
#include "Iterator.h"
#include "IOUtil.h"
#include "MemoryUtil.h"
#include "SAM.h"
#include "StringUtil.h" // for toSI
#include "Uncompress.h"
#include "Pipe.h"
#include "PipeMux.h"
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
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/time.h>

using namespace std;

#define PROGRAM "KAligner"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2012 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... QUERY... TARGET\n"
"Align the sequences of the files QUERY to those of TARGET.\n"
"All perfect matches of at least k bases will be found.\n"
"\n"
"  -k, -l, --kmer=N      k-mer size and minimum alignment length\n"
"  -s, --section=S/N     split the target into N sections and align"
"                        reads to section S [1/1]\n"
"  -i, --ignore-multimap ignore duplicate k-mer in the target\n"
"                        [default]\n"
"  -m, --multimap        allow duplicate k-mer in the target\n"
"      --no-multimap     disallow duplicate k-mer in the target\n"
"  -j, --threads=N       use N threads [2] up to one per query file\n"
"                        or if N is 0 use one thread per query file\n"
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
	static unsigned section = 1;
	static unsigned nsections = 1;

	/** Output formats */
	static int format;
}

static const char shortopts[] = "ij:k:l:mo:s:v";


enum { OPT_HELP = 1, OPT_VERSION, OPT_SYNC };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "section",     required_argument, NULL, 's' },
	{ "no-multi",    no_argument,     &opt::multimap, opt::ERROR },
	{ "multimap",    no_argument,     &opt::multimap, opt::MULTIMAP },
	{ "ignore-multimap", no_argument, &opt::multimap, opt::IGNORE },
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
				cerr << "error: unexpected character: "
					"`" << c << "'\n";
				exit(EXIT_FAILURE);
			}
			break;
		}
	}

	size_t overlaps = contigs * (opt::k-1);
	size_t kmer = bases - overlaps;
	if (opt::verbose > 0) {
		cerr << "Read " << bases << " bases, "
			<< contigs << " contigs, " << scaffolds << " scaffolds"
			" from `" << path << "'. "
			"Expecting " << kmer << " k-mer.\n";
		cerr << "Index will use at least "
			<< toSI(kmer * sizeof(pair<Kmer, Position>))
			<< "B.\n";
	}
	assert(bases > overlaps);
	return kmer;
}

template <class SeqPosHashMap>
static void readContigsIntoDB(string refFastaFile,
		Aligner<SeqPosHashMap>& aligner);
static void *alignReadsToDB(void *arg);
static void *readFile(void *arg);

/** Unique aligner using map */
static Aligner<SeqPosHashUniqueMap> *g_aligner_u;

/** Multimap aligner using multimap */
static Aligner<SeqPosHashMultiMap> *g_aligner_m;

/** Number of reads. */
static unsigned g_readCount;

/** Number of reads that aligned. */
static unsigned g_alignedCount;

/** Guard cerr. */
static pthread_mutex_t g_mutexCerr = PTHREAD_MUTEX_INITIALIZER;

/** Stores the output string and the read index number for an
 * alignment. */
struct OutData
{
	string s;
	size_t index;

	OutData(string s = string(), size_t index = 0)
		: s(s), index(index) { }

	/** Operator needed for sorting priority queue. */
	bool operator<(const OutData& a) const
	{
		// Smaller index number has higher priority.
		return index > a.index;
	}
};

/** Shares data between workers and the output thread. */
static Pipe<OutData> g_pipeOut(1<<7);

/** Shares data between producer and worker threads. */
static PipeMux<FastaRecord> g_pipeMux(1);

/** Notification of the current size of the g_pqueue. */
static size_t g_pqSize;
static const size_t MAX_PQ_SIZE = 1000;

/** Conditional variable used to block workers until the g_pqueue has
 * become small enough. */
static pthread_cond_t g_pqSize_cv = PTHREAD_COND_INITIALIZER;
static pthread_mutex_t g_mutexPqSize = PTHREAD_MUTEX_INITIALIZER;

static void* printAlignments(void*)
{
	priority_queue<OutData> pqueue;
	size_t index = 1;
	for (pair<OutData, size_t> p = g_pipeOut.pop();
			p.second > 0; p = g_pipeOut.pop()) {
		pqueue.push(p.first);
		while (!pqueue.empty()) {
			const OutData& rec = pqueue.top();
			if (index == rec.index) {
				// Print the record at the current index.
				index++;
				assert(rec.index > 0);
				cout << rec.s;
				assert_good(cout, "stdout");
				pqueue.pop();
			} else if (g_pipeMux.invalidEntry(index)) {
				// Skip this index since it is invalid.
				index++;
			} else {
				// The record for this index has not been added, get
				// another record from the pipe.
				break;
			}
		}

		// Let waiting workers continue if the pqueue is small enough.
		pthread_mutex_lock(&g_mutexPqSize);
		g_pqSize = pqueue.size();
		if (g_pqSize < MAX_PQ_SIZE)
			pthread_cond_broadcast(&g_pqSize_cv);
		pthread_mutex_unlock(&g_mutexPqSize);
	}
	return NULL;
}

/** Store a FastaReader and Pipe. */
struct WorkerArg {
	FastaReader& in;
	Pipe<FastaRecord>& out;
	WorkerArg(FastaReader& in, Pipe<FastaRecord>& out)
		: in(in), out(out) { }
	~WorkerArg() { delete &in; }
};

static pthread_t getReadFiles(const char *readsFile)
{
	if (opt::verbose > 0) {
		pthread_mutex_lock(&g_mutexCerr);
		cerr << "Reading `" << readsFile << "'...\n";
		pthread_mutex_unlock(&g_mutexCerr);
	}

	FastaReader* in = new FastaReader(
			readsFile, FastaReader::FOLD_CASE);
	WorkerArg* arg = new WorkerArg(*in, *g_pipeMux.addPipe());

	pthread_t thread;
	pthread_create(&thread, NULL, readFile, static_cast<void*>(arg));

	return thread;
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

	char delim = '/';
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'k': case 'l':
				arg >> opt::k;
				break;
			case 'm': opt::multimap = opt::MULTIMAP; break;
			case 'i': opt::multimap = opt::IGNORE; break;
			case 'j': arg >> opt::threads; break;
			case 'v': opt::verbose++; break;
			case 's': arg >> opt::section >> delim >>
					  opt::nsections; break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
		if (optarg != NULL && !arg.eof()) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	if (opt::section == 0 || opt::section > opt::nsections
			|| delim != '/') {
		cerr << PROGRAM ": -s, --section option is incorrectly set\n";
		die = true;
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

	g_readCount = 0;

	vector<pthread_t> producer_threads;
	transform(argv + optind, argv + argc,
			back_inserter(producer_threads), getReadFiles);

	vector<pthread_t> threads;
	for (int i = 0; i < opt::threads; i++) {
		pthread_t thread;
		pthread_create(&thread, NULL, alignReadsToDB, NULL);
		threads.push_back(thread);
	}

	pthread_t out_thread;
	pthread_create(&out_thread, NULL, printAlignments, NULL);

	void *status;
	// Wait for all threads to finish.
	for (size_t i = 0; i < producer_threads.size(); i++)
		pthread_join(producer_threads[i], &status);
	for (size_t i = 0; i < threads.size(); i++)
		pthread_join(threads[i], &status);
	g_pipeOut.close();
	pthread_join(out_thread, &status);

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
		"Hash load: " << size << " / " << buckets
		<< " = " << (float)size / buckets
		<< " using " << toSI(getMemoryUsage()) << "B." << endl;
}

template <class SeqPosHashMap>
static void readContigsIntoDB(string refFastaFile,
		Aligner<SeqPosHashMap>& aligner)
{
	if (opt::verbose > 0)
		cerr << "Reading target `" << refFastaFile << "'..." << endl;

	unsigned count = 0;
	FastaReader in(refFastaFile.c_str(), FastaReader::FOLD_CASE);
	if (opt::nsections > 1)
		in.split(opt::section, opt::nsections);
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

		cout << "@SQ\tSN:" << rec.id
			<< "\tLN:" << rec.seq.length() << '\n';
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

/** Read each fasta record from 'in', and add it to 'pipe'. */
static void readFile(FastaReader& in, Pipe<FastaRecord>& pipe)
{
	for (FastaRecord rec; in >> rec;)
		pipe.push(rec);
	assert(in.eof());
	pipe.close();
}

/** Producer thread. */
static void* readFile(void* arg)
{
	WorkerArg* p = static_cast<WorkerArg*>(arg);
	readFile(p->in, p->out);
	delete p;
	return NULL;
}

/** @Returns the time in seconds between [start, end]. */
static double timeDiff(const timeval& start, const timeval& end)
{
	double result = (double)end.tv_sec +
		(double)end.tv_usec/1000000.0;
	result -= (double)start.tv_sec +
		(double)start.tv_usec/1000000.0;
	return result;
}

static void* alignReadsToDB(void*)
{
	opt::chastityFilter = false;
	opt::trimMasked = false;
	static timeval start, end;

	pthread_mutex_lock(&g_mutexCerr);
	gettimeofday(&start, NULL);
	pthread_mutex_unlock(&g_mutexCerr);

	for (pair<FastaRecord, size_t> recPair = g_pipeMux.nextValue();
			recPair.second > 0; recPair = g_pipeMux.nextValue()) {
		const FastaRecord& rec = recPair.first;
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

		ostringstream out;
		string s = output.str();
		switch (opt::format) {
		  case KALIGNER:
			out << rec.id;
			if (opt::printSeq) {
				out << ' ';
				if (opt::colourSpace)
					out << rec.anchor;
				out << seq;
			}
			out << s << '\n';
			break;
		  case SAM:
			out << s;
			break;
		}
		g_pipeOut.push(OutData(out.str(), recPair.second));

		// Prevent the priority_queue from growing too large by
		// waiting for threads going far too slow.
		pthread_mutex_lock(&g_mutexPqSize);
		if (g_pqSize >= MAX_PQ_SIZE)
			pthread_cond_wait(&g_pqSize_cv, &g_mutexPqSize);
		pthread_mutex_unlock(&g_mutexPqSize);

		if (opt::verbose > 0) {
			pthread_mutex_lock(&g_mutexCerr);
			if (!s.empty())
				g_alignedCount++;
			if (++g_readCount % 1000000 == 0) {
				gettimeofday(&end, NULL);
				double result = timeDiff(start, end);
				cerr << "Aligned " << g_readCount << " reads at "
					<< (int)(1000000 / result) << " reads/sec.\n";
				start = end;
			}
			pthread_mutex_unlock(&g_mutexCerr);
		}
	}
	return NULL;
}
