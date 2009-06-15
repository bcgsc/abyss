#include "Aligner.h"
#include "PairedAlgorithms.h"
#include "FastaReader.h"
#include <cassert>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <list>
#include <sstream>
#include <string>

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
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static unsigned k;
	static int verbose;
}

static const char* shortopts = "k:o:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

static void readContigsIntoDB(string refFastaFile, Aligner& aligner);
static void alignReadsToDB(string readsFile, Aligner& aligner);

int main(int argc, char** argv)
{
	bool die = false;
	for (char c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
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
	} else if (argc - optind > 2) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	string readsFile(argv[optind++]);
	string refFastaFile(argv[optind++]);

	if (opt::verbose > 0)
		cerr << "k: " << opt::k
			<< " Query: " << readsFile
			<< " Target: " << refFastaFile
			<< endl;

	Aligner aligner(opt::k);
	readContigsIntoDB(refFastaFile, aligner);
	alignReadsToDB(readsFile, aligner);
	return 0;
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
		aligner.addReferenceSequence(contigID, seq);

		count++;
		if (opt::verbose > 0 && count % 100000 == 0)
			printProgress(aligner, count);
	}
	if (opt::verbose > 0)
			printProgress(aligner, count);

	fileHandle.close();
}

static void alignReadsToDB(string readsFile, Aligner& aligner)
{
	FastaReader fileHandle(readsFile.c_str());

	unsigned count = 0;
	while (fileHandle.isGood()) {
		string readID;
		Sequence readSeq = fileHandle.ReadSequence(readID);

		list<Alignment> alignments;
		size_t pos = readSeq.find_first_not_of("ACGT");
		if (pos == string::npos)
			aligner.alignRead(readSeq, back_inserter(alignments));

		cout << readID;
		for (list<Alignment>::iterator iter = alignments.begin();
				iter != alignments.end(); ++iter)
			cout << '\t' << *iter;
		cout << '\n';
		assert(cout.good());

		if (opt::verbose > 0 && ++count % 1000000 == 0)
			cerr << "Aligned " << count << " reads\n";
	}
	if (opt::verbose > 0)
		cerr << "Aligned " << count << " reads\n";
}
