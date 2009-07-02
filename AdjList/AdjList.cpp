#include "FastaReader.h"
#include "PackedSeq.h"
#include "PairUtils.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring> // for strerror
#include <fstream>
#include <functional>
#include <iterator>
#include <getopt.h>
#include <iostream>

using namespace std;

#define PROGRAM "AdjList"

static const char *VERSION_MESSAGE =
PROGRAM " (ABySS) " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2009 Canada's Michael Smith Genome Science Centre\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION]... [FILE]...\n"
"Find all contigs that overlap by exactly k-1 bases. Contigs may be read\n"
"from FILE(s) or standard input. Output is written to standard output.\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static int k;
	static int verbose;
	extern bool colourSpace;
}

static const char* shortopts = "k:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",    required_argument, NULL, 'k' },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

void generatePossibleExtensions(const PackedSeq& seq, extDirection dir, PSequenceVector& outseqs);

/** A contig ID and its end sequences. */
struct ContigEndSeq {
	ContigID id;
	PackedSeq l;
	PackedSeq r;
	ContigEndSeq(const ContigID& id,
			const PackedSeq& l, const PackedSeq& r)
		: id(id), l(l), r(r) { }
};

static void readContigs(string path, vector<ContigEndSeq>* pContigs)
{
	unsigned count = 0;
	FastaReader in(path.c_str());
	while (in.isGood()) {
		ContigID id;
		Sequence seq = in.ReadSequence(id);

		if (count++ == 0) {
			// Detect colour-space contigs.
			opt::colourSpace = isdigit(seq[0]);
		} else {
			if (opt::colourSpace)
				assert(isdigit(seq[0]));
			else
				assert(isalpha(seq[0]));
		}

		PackedSeq seql = seq.substr(seq.length() - opt::k, opt::k);
		PackedSeq seqr = seq.substr(0, opt::k);
		pContigs->push_back(ContigEndSeq(id, seql, seqr));
	}
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

	vector<ContigEndSeq> contigs;
	if (optind < argc) {
		for_each(argv + optind, argv + argc,
				bind2nd(ptr_fun(readContigs), &contigs));
	} else
		readContigs("-", &contigs);

	// Generate a k-mer -> contig lookup table for all the contig ends
	std::map<PackedSeq, ContigID> contigLUTs[2];
	for (vector<ContigEndSeq>::const_iterator i = contigs.begin();
			i != contigs.end(); ++i) {
		contigLUTs[0][i->l] = i->id;
		contigLUTs[1][i->r] = i->id;
	}

	ostream& out = cout;
	int numVerts = 0;
	int numEdges = 0;
	for (vector<ContigEndSeq>::const_iterator i = contigs.begin();
			i != contigs.end(); ++i) {
		const ContigID& id = i->id;
		const PackedSeq seqs[2] = { i->l, i->r };

		out << id;

		const unsigned numEnds = 2;
		for(unsigned idx = 0; idx < numEnds; idx++)
		{
			std::vector<SimpleEdgeDesc> edges;
			const PackedSeq& currSeq = seqs[idx];
			extDirection dir;
			dir = (idx == 0) ? SENSE : ANTISENSE;

			
			// Generate the links
			PSequenceVector extensions;
			generatePossibleExtensions(currSeq, dir, extensions);
		  
			for(PSequenceVector::iterator iter = extensions.begin(); iter != extensions.end(); ++iter)
			{
				// Get the contig this sequence maps to
				for(size_t compIdx = 0; compIdx <= 1; ++compIdx)
				{
					size_t lookuptable_id = 1 - idx;
					bool reverse = (compIdx == 1);
					PackedSeq testSeq;
					if(reverse)
					{
						// flip the lookup table id
						lookuptable_id = 1 - lookuptable_id;
						testSeq = reverseComplement(*iter);
					}
					else
					{
						testSeq = *iter;
					}
					std::map<PackedSeq, ContigID>::iterator cLUTIter;
					cLUTIter = contigLUTs[lookuptable_id].find(testSeq);
					if(cLUTIter != contigLUTs[lookuptable_id].end())
					{						
						// Store the edge
						SimpleEdgeDesc ed;
						ed.contig = cLUTIter->second;
						ed.isRC = reverse;
						edges.push_back(ed);
					}
				}
			}
			// Print the edges
			out << " [ ";
			std::copy(edges.begin(), edges.end(), std::ostream_iterator<SimpleEdgeDesc>(out, " "));
			out << ']';
			numEdges += edges.size();
		}
		out << '\n';
		numVerts++;
	}

	if (opt::verbose > 0)
		cerr << "vertices: " << numVerts << " "
			"edges: " << numEdges << endl;
} 

void generatePossibleExtensions(const PackedSeq& seq, extDirection dir, PSequenceVector& outseqs)
{
  PackedSeq extSeq(seq);
  extSeq.shift(dir);
  
  // Check for the existance of the 4 possible extensions
  for(int i  = 0; i < NUM_BASES; i++)
    {
      extSeq.setLastBase(dir, i);
      outseqs.push_back(extSeq);
    }
}
