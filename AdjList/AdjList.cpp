#include "Common/Options.h"
#include "FastaReader.h"
#include "HashMap.h"
#include "Kmer.h"
#include "Uncompress.h"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <cstring> // for strerror
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

using namespace std;

#define PROGRAM "AdjList"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2009 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [FILE]...\n"
"Find all contigs that overlap by exactly k-1 bases. Contigs may be read\n"
"from FILE(s) or standard input. Output is written to standard output.\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"      --adj             output the results in adj format [DEFAULT]\n"
"      --dot             output the results in dot format\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

/** Enumeration of output formats */
enum format { ADJ, DOT };

namespace opt {
	static int k;
	static int overlap;

	/** Output formats */
	static int format;
}

static const char shortopts[] = "k:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",    required_argument, NULL, 'k' },
	{ "adj",     no_argument,       &opt::format, ADJ },
	{ "dot",     no_argument,       &opt::format, DOT },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

typedef string ContigID;

/** A contig ID and its end sequences. */
struct ContigEndSeq {
	ContigID id;
	unsigned length;
	Kmer l;
	Kmer r;
	ContigEndSeq(const ContigID& id, unsigned length,
			const Kmer& l, const Kmer& r)
		: id(id), length(length), l(l), r(r) { }
};

static void readContigs(string path, vector<ContigEndSeq>* pContigs)
{
	if (opt::verbose > 0)
		cerr << "Reading `" << path << "'...\n";

	unsigned count = 0;
	FastaReader in(path.c_str(), FastaReader::KEEP_N);
	for (FastaRecord rec; in >> rec;) {
		const Sequence& seq = rec.seq;
		if (count++ == 0) {
			// Detect colour-space contigs.
			opt::colourSpace = isdigit(seq[0]);
		} else {
			if (opt::colourSpace)
				assert(isdigit(seq[0]));
			else
				assert(isalpha(seq[0]));
		}

		Kmer seql(seq.substr(seq.length() - opt::overlap,
				opt::overlap));
		Kmer seqr(seq.substr(0, opt::overlap));
		pContigs->push_back(ContigEndSeq(rec.id, seq.length(),
					seql, seqr));
	}
	assert(in.eof());
}

/** The relevant information of our contigs.
 * This variable is global so that ContigNode may use the mapping of
 * contig numbers to contig names. */
static vector<ContigEndSeq> contigs;

struct ContigNode {
	unsigned id;
	bool sense;

	ContigNode(unsigned id, bool sense) : id(id), sense(sense) { }

	friend ostream& operator <<(ostream& out, const ContigNode& o)
	{
		return out << contigs[o.id].id << ',' << o.sense;
	}
};

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
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
		cerr << PROGRAM ": " << "missing -k,--kmer option\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	opt::overlap = opt::k - 1;
	Kmer::setLength(opt::overlap);

	if (optind < argc) {
		for_each(argv + optind, argv + argc,
				bind2nd(ptr_fun(readContigs), &contigs));
	} else
		readContigs("-", &contigs);

	if (opt::verbose > 0)
		cerr << "Read " << contigs.size() << " contigs\n";

	typedef hash_map<Kmer, vector<ContigNode>, hashKmer> KmerMap;
	vector<KmerMap> ends(2, KmerMap(contigs.size()));
	for (vector<ContigEndSeq>::const_iterator it = contigs.begin();
			it != contigs.end(); ++it) {
		unsigned i = it - contigs.begin();
		ends[0][it->l].push_back(ContigNode(i, false));
		ends[1][reverseComplement(it->l)].push_back(
				ContigNode(i, true));
		ends[1][it->r].push_back(ContigNode(i, false));
		ends[0][reverseComplement(it->r)].push_back(
				ContigNode(i, true));
	}

	ostream& out = cout;
	if (opt::format == DOT)
		out << "digraph adj {\n";

	int numVerts = 0;
	int numEdges = 0;
	for (vector<ContigEndSeq>::const_iterator i = contigs.begin();
			i != contigs.end(); ++i) {
		const ContigID& id = i->id;

		if (opt::format == ADJ)
			out << id << ' ' << i->length;

		for (unsigned idx = 0; idx < 2; idx++) {
			const Kmer& seq = idx == 0 ? i->l : i->r;
			const KmerMap::mapped_type& edges = ends[!idx][seq];

			switch (opt::format) {
			  case ADJ:
				out << " [ ";
				copy(edges.begin(), edges.end(),
						ostream_iterator<ContigNode>(out, " "));
				out << ']';
				break;
			  case DOT:
				out << '"' << id << (idx ? '-' : '+') << "\" [len="
					<< i->length << "];\n"
					<< '"' << id << (idx ? '-' : '+') << '"';
				if (!edges.empty()) {
					out << " -> {";
					for (KmerMap::mapped_type::const_iterator it
							= edges.begin(); it != edges.end(); ++it)
						out << " \"" << contigs[it->id].id
							<< (idx != it->sense ? '-' : '+') << '"';
					out << " }";
				}
				out << ";\n";
				break;
			}
			numEdges += edges.size();
		}
		if (opt::format == ADJ)
			out << '\n';
		numVerts++;
	}

	if (opt::format == DOT)
		out << "}\n";

	if (opt::verbose > 0)
		cerr << "Vertices: " << numVerts
			<< " Edges: " << numEdges << endl;
}
