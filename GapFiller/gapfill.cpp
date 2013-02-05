#include "SAM.h"
#include "FastaReader.h"

#include "config.h"
#include <cstdlib>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <map>


using namespace std;

#define PROGRAM "abyss-gapfill"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Anthony Raymond.\n"
"\n"
"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";

//TODO
static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... ALIGNS CONTIGS\n"
"Write read pairs that map to the same contig to the file SAME.\n"
"Write read pairs that map to different contigs to stdout.\n"
"Alignments may be in ALIGNS or standard input.\n"
"\n"
"      --no-qname        set the qname to * [default]\n"
"      --qname           do not alter the qname\n"
"  -l, --min-align=N     the minimal alignment size [1]\n"
"  -s, --same=SAME       write properly-paired reads to this file\n"
"  -h, --hist=FILE       write the fragment size histogram to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

//TODO
namespace opt {
	static string alignPath;
	static int verbose;
}

//TODO
static const char shortopts[] = "h:l:s:v";

//TODO
enum { OPT_HELP = 1, OPT_VERSION };

//TODO
static const struct option longopts[] = {
	{ "min-align", required_argument, NULL, 'l' },
	{ "hist",    required_argument, NULL, 'h' },
	{ "same",    required_argument, NULL, 's' },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

typedef map<string, FastaRecord> Scaffolds;
typedef multimap<string, SAMRecord> Alignments;

static struct {
	int seqs;
	int scaffolds;
	int aligns;
	int split_same;
} stats;

static void readScaffolds(string path, Scaffolds& scaffs)
{
	FastaReader in(path.c_str(), FastaReader::NO_FOLD_CASE);
	FastaRecord rec;
	while (in >> rec) {
		// Store only scaffolded sequences
		stats.seqs++;
		size_t i = rec.seq.find_first_of("N");
		if (i != string::npos) {
			stats.scaffolds++;
			assert(scaffs.find(rec.id) == scaffs.end());
			scaffs.insert(make_pair(rec.id, rec));
		}
	}
}

static void readAlignments(string path, Alignments& aligns)
//		Scaffolds& scaffs)
{
	ifstream in(path.c_str());
	SAMRecord rec;
	while (in.peek() == '@') {
		string tmp;
		getline(in, tmp);
		assert(in);
	}

	while (in >> rec) {
		// If aligns to either side of gap store in aligns.
		// Lets start with just reads that align to the same contig...
		stats.aligns++;
		if (rec.tags.size() >= 2 && rec.tags[0] == 'X'
				&& rec.tags[1] == 'A'
				&& rec.tags.find(rec.rname, 2) != string::npos) {
			stats.split_same++;
			aligns.insert(make_pair(rec.rname, rec));
		}
		assert(in);
	}
	assert(in.eof());
}

int main(int argc, char* const* argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'v': opt::verbose++; break;
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

	if (argc - optind != 2) {
		cerr << PROGRAM ": incorrect number of arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	Scaffolds scaffs;
	readScaffolds(argv[optind++], scaffs);
	cerr << "sequences: " << stats.seqs << "\tscaffolds: "
		<< stats.scaffolds << '\n';

	Alignments aligns;
	readAlignments(argv[optind++], aligns);
	cerr << "aligns: " << stats.aligns << "\tsplit: "
		<< stats.split_same << '\n';

	

	//fillGaps();
	//writeContigs();

	return 0;
}
