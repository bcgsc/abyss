#include "config.h"
#include "Dictionary.h"
#include "FastaReader.h"
#include <algorithm>
#include <cstdlib>
#include <cerrno>
#include <cstring> // for strerror
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#define PROGRAM "MergeContigs"

static const char *VERSION_MESSAGE =
PROGRAM " (ABySS) " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2009 Canada's Michael Smith Genome Science Centre\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION]... CONTIG PATH\n"
"Merge paths of contigs to create larger contigs.\n"
"  CONTIG  contigs in FASTA format\n"
"  PATH    paths of these contigs\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -o, --out=FILE        write result to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static unsigned k;
	static int verbose;
	static string out;
	extern bool colourSpace;
}

static const char* shortopts = "k:o:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "out",         required_argument, NULL, 'o' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Return the last character of s and remove it. */
static char chop(string& s)
{
	assert(s.length() > 1);
	unsigned back = s.length() - 1;
	char c = s[back];
	s.erase(back);
	return c;
}

static unsigned line_num;

static void assert_plus_minus(char c)
{
	if (c != '+' && c != '-') {
		cerr << "error: " << line_num
			<< ": expected `+' or `-' and saw `" << c << '\''
			<< endl;
		exit(EXIT_FAILURE);
	}
}

struct ContigNode {
	string id;
	bool sense;

	friend istream& operator >>(istream& in, ContigNode& o)
	{
		if (in >> o.id) {
			char c = chop(o.id);
			assert_plus_minus(c);
			o.sense = c == '-';
		}
		return in;
	}

	friend ostream& operator <<(ostream& out, const ContigNode& o)
	{
		return out << o.id << (o.sense ? '-' : '+');
	}
};

struct Path : vector<ContigNode>
{
	friend ostream& operator <<(ostream& out, const Path& o)
	{
		copy(o.begin(), o.end()-1,
				ostream_iterator<ContigNode>(out, ","));
		return out << o.back();
	}
};

static const string& contigNodeToId(const ContigNode& node)
{
	return node.id;
}

static set<string> getContigIDs(const vector<Path>& paths)
{
	set<string> seen;
	for (vector<Path>::const_iterator it = paths.begin();
			it != paths.end(); ++it)
		transform(it->begin(), it->end(),
				inserter(seen, seen.begin()),
				contigNodeToId);
	return seen;
}

static void assert_open(std::ifstream& f, const std::string& p)
{
	if (f.is_open())
		return;
	std::cerr << p << ": " << strerror(errno) << std::endl;
	exit(EXIT_FAILURE);
}

static Dictionary g_dict;

struct Contig {
	string id;
    Sequence seq;
    unsigned coverage;
    Contig(string id, Sequence seq, unsigned coverage)
        : id(id), seq(seq), coverage(coverage) { }

	operator FastaRecord()
	{
		ostringstream s;
		s << seq.length() << ' ' << coverage;
		return FastaRecord(id, s.str(), seq);
	}

	friend ostream& operator <<(ostream& out, const Contig& o)
	{
		return out << '>' << o.id << ' '
			<< o.seq.length() << ' ' << o.coverage << '\n'
			<< o.seq << '\n';
	}
};

static Contig mergePath(const Path& path,
		const vector<Contig>& contigs)
{
	Sequence seq;
	unsigned coverage = 0;
	for (Path::const_iterator it = path.begin();
			it != path.end(); ++it) {
		unsigned serial = g_dict.serial(it->id);
		const Contig& contig = contigs[serial];
		coverage += contig.coverage;

		Sequence h = it->sense ? reverseComplement(contig.seq)
			: contig.seq;
		if (seq.empty()) {
			seq = h;
			continue;
		}

		unsigned overlap = opt::k - 1;
		Sequence a(seq, 0, seq.length() - overlap);
		Sequence ao(seq, seq.length() - overlap);
		Sequence bo(h, 0, overlap);
		Sequence b(h, overlap);
		if (ao != bo) {
			cerr << "error: the head of `" << contig.id << "' "
				"does not match the tail of the previous contig\n"
				<< ao << '\n' << bo << '\n' << path << endl;
			exit(EXIT_FAILURE);
		}
		seq += b;
	}
	return Contig("", seq, coverage);
}

template<typename T> static string toString(T x)
{
	ostringstream s;
	s << x;
	return s.str();
}

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
			case 'o': arg >> opt::out; break;
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

	if (opt::out.empty()) {
		cerr << PROGRAM ": " << "missing -o,--out option\n";
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

	const char* contigFile = argv[optind++];
	string pathFile(argv[optind++]);

	vector<Contig> contigs;
	{
		FastaReader in(contigFile, FastaReader::KEEP_N);
		for (FastaRecord rec; in >> rec;) {
			istringstream ss(rec.comment);
			unsigned length, coverage = 0;
			ss >> length >> coverage;
			unsigned serial = g_dict.serial(rec.id);
			assert(contigs.size() == serial);
			contigs.push_back(Contig(rec.id, rec.seq, coverage));
		}
		assert(in.eof());
		assert(!contigs.empty());
		opt::colourSpace = isdigit(contigs[0].seq[0]);
	}

	vector<Path> paths;
	{
		ifstream fin(pathFile.c_str());
		if (pathFile != "-")
			assert_open(fin, pathFile);
		istream& in = pathFile == "-" ? cin : fin;
		for (string s; getline(in, s);) {
			line_num++;
			istringstream ss(s);
			Path path;
			copy(istream_iterator<ContigNode>(ss),
					istream_iterator<ContigNode>(),
					back_inserter(path));
			paths.push_back(path);
		}
		assert(in.eof());
	}

	ofstream out(opt::out.c_str());
	assert(out.good());

	set<string> seen = getContigIDs(paths);
	for (vector<Contig>::const_iterator it = contigs.begin();
			it != contigs.end(); ++it)
		if (seen.count(it->id) == 0)
			out << *it;

	int id = contigs.size();
	for (vector<Path>::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		const Path& path = *it;
		Contig contig = mergePath(path, contigs);
		contig.id = toString(id++);
		FastaRecord rec(contig);
		rec.comment += ' ' + toString(path);
		out << rec;
		assert(out.good());
	}

	return 0;
}
