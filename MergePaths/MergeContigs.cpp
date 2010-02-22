#include "config.h"
#include "Common/Options.h"
#include "ContigNode.h"
#include "Dictionary.h"
#include "FastaReader.h"
#include "StringUtil.h"
#include <algorithm>
#include <cstdlib>
#include <cerrno>
#include <cstring> // for strerror
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#define PROGRAM "MergeContigs"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2010 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... CONTIG PATH\n"
"Merge paths of contigs to create larger contigs.\n"
"  CONTIG  contigs in FASTA format\n"
"  PATH    paths of these contigs\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -o, --out=FILE        write result to FILE\n"
"  -p, --path=PATH_FILE  paths output by SimpleGraph\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static unsigned k;
	static string out;
	static string path;
}

static const char shortopts[] = "k:o:p:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "out",         required_argument, NULL, 'o' },
	{ "path",        required_argument, NULL, 'p' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
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

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

struct Contig {
	string id;
    Sequence seq;
    unsigned coverage;
    Contig(const string& id, const Sequence& seq, unsigned coverage)
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

static vector<Contig> g_contigs;

/** Return the sequence of the specified contig node. The sequence
 * may be ambiguous or reverse complemented.
 */
const Sequence ContigNode::sequence() const
{
	if (ambiguous()) {
		return string(opt::k - 1, 'N') + ambiguousSequence();
	} else {
		const Sequence& seq = g_contigs[id()].seq;
		return sense() ? reverseComplement(seq) : seq;
	}
}

/** Return a consensus sequence of a and b.
 * @return an empty string if a consensus could not be found
 */
static string createConsensus(const Sequence& a, const Sequence& b)
{
	assert(a.length() == b.length());
	string s;
	s.reserve(a.length());
	for (string::const_iterator ita = a.begin(), itb = b.begin();
			ita != a.end(); ++ita, ++itb) {
		char c = *ita == *itb ? *ita
			: *ita == 'N' ? *itb
			: *itb == 'N' ? *ita
			: 'x';
		if (c == 'x')
			return string("");
		s += c;
	}
	return s;
}

/** Merge the specified two contigs, which must overlap by k-1 bp, and
 * generate a consensus sequence of the overlapping region. The result
 * is stored in the first argument.
 */
static void mergeContigs(Sequence& seq, const Sequence& s,
		const ContigNode& node, const Path& path)
{
	unsigned overlap = opt::k - 1;
	assert(seq.length() > overlap);
	assert(s.length() > overlap);
	Sequence ao(seq, seq.length() - overlap);
	Sequence bo(s, 0, overlap);
	Sequence o = createConsensus(ao, bo);
	if (o.empty()) {
		cerr << "error: the head of `" << node << "' "
			"does not match the tail of the previous contig\n"
			<< ao << '\n' << bo << '\n' << path << endl;
		exit(EXIT_FAILURE);
	}
	seq.resize(seq.length() - overlap);
	seq += o;
	seq += Sequence(s, overlap);
}

static Contig mergePath(const Path& path,
		const vector<Contig>& contigs)
{
	Sequence seq;
	unsigned coverage = 0;
	for (Path::const_iterator it = path.begin();
			it != path.end(); ++it) {
		if (!it->ambiguous())
			coverage += contigs[it->id()].coverage;
		if (seq.empty())
			seq = it->sequence();
		else
			mergeContigs(seq, it->sequence(), *it, path);
	}
	return Contig("", seq, coverage);
}

template<typename T> static string toString(T x)
{
	ostringstream s;
	s << x;
	return s.str();
}

/** Loads all paths from the file named inPath into paths. */
static void loadPaths(string& inPath, vector<Path>& paths)
{
	ifstream fin(inPath.c_str());
	if (opt::verbose > 0)
		cerr << "Reading `" << inPath << "'..." << endl;
	if (inPath != "-")
		assert_open(fin, inPath);
	istream& in = inPath == "-" ? cin : fin;

	for (string s; getline(in, s);) {
		istringstream ss(s);

		string sID;
		ss >> sID;

		Path path;
		copy(istream_iterator<ContigNode>(ss),
				istream_iterator<ContigNode>(),
				back_inserter(path));
		paths.push_back(path);
	}
	assert(in.eof());
}

/** Finds all contigs used in each path in paths, and
 * marks them as seen in the vector seen. */
static void seenContigs(vector<bool>& seen, const vector<Path>& paths)
{
	for (vector<Path>::const_iterator it = paths.begin();
			it != paths.end(); ++it)
		for (Path::const_iterator itc = it->begin();
				itc != it->end(); ++itc)
			if (itc->id() < seen.size())
				seen[itc->id()] = true;
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
			case 'p': arg >> opt::path; break;
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
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	const char* contigFile = argv[optind++];
	string mergedPathFile(argv[optind++]);

	vector<Contig>& contigs = g_contigs;
	{
		FastaReader in(contigFile, FastaReader::KEEP_N);
		for (FastaRecord rec; in >> rec;) {
			istringstream ss(rec.comment);
			unsigned length, coverage = 0;
			ss >> length >> coverage;
			unsigned serial = g_contigIDs.serial(rec.id);
			assert(contigs.size() == serial);
			(void)serial;
			contigs.push_back(Contig(rec.id, rec.seq, coverage));
		}
		assert(in.eof());
		assert(!contigs.empty());
		opt::colourSpace = isdigit(contigs[0].seq[0]);
		if (argc - optind == 0) g_contigIDs.lock();
	}

	vector<Path> paths;
	loadPaths(mergedPathFile, paths); 
	if (opt::verbose > 0)
		cerr << "Total number of paths: " << paths.size() << '\n';

	// Record all the contigs that are in a path.
	vector<bool> seen(contigs.size());
	seenContigs(seen, paths);

	// Record all the contigs that were in a previous path.
	if (argc - optind > 0) {
		unsigned count = 0;
		for (; optind < argc; optind++) {
			vector<Path> prevPaths;
			string filename(argv[optind]);
			loadPaths(filename, prevPaths);
			seenContigs(seen, prevPaths);
			count += prevPaths.size();
		}
		if (opt::verbose > 0)
			cerr << "Total number of previous paths: " << count << '\n';
	}

	// Record all the contigs that are seeds.
	if (!opt::path.empty()) {
		vector<bool> seenPivots(contigs.size());
		ifstream fin(opt::path.c_str());
		assert_open(fin, opt::path);
		string s;
		while (fin >> s) {
			fin.ignore(numeric_limits<streamsize>::max(), '\n');
			unsigned pivotNum = g_contigIDs.serial(s);
			assert(pivotNum < contigs.size());
			// Only count a pivot as seen if it was in a final path.
			if (seen[pivotNum])
				seenPivots[pivotNum] = true;
		}
		assert(fin.eof());
		seen = seenPivots;
	}

	// Output those contigs that were not seen in a path.
	ofstream out(opt::out.c_str());
	assert(out.good());
	for (vector<Contig>::const_iterator it = contigs.begin();
			it != contigs.end(); ++it)
		if (!seen[g_contigIDs.serial(it->id)])
			out << *it;

	int id;
	stringstream s(g_contigIDs.key(contigs.size() - 1));
	s >> id;
	id++;
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
