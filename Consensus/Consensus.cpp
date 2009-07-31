#include "Aligner.h"
#include "FastaReader.h"
#include "PairUtils.h"
#include "FastaWriter.h"
#include <getopt.h>
#include <string>
#include <sstream>
#include <iostream>

using namespace std;

#define PROGRAM "ParseAligns"

static const char *VERSION_MESSAGE =
PROGRAM " (ABySS) " VERSION "\n"
"Written by Shaun Jackman and Tony Raymond.\n"
"\n"
"Copyright 2009 Canada's Michael Smith Genome Science Centre\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION]... [FILE]...\n"
"Write read pairs that align to the same contig to FRAGMENTS or HISTOGRAM.\n"
"Write read pairs that align to different contigs to standard output.\n"
"Alignments may be in FILE(s) or standard input.\n"
"\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static string outPath;
	static int verbose;
	extern bool colourSpace;
	static bool csToNt;
}

static const char* shortopts = "v:o:C";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "out",		 required_argument, NULL, 'o' },
	{ "cstont",		 no_argument,		NULL, 'C' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

struct BaseCount {
	unsigned count[4];
	friend ostream& operator <<(ostream& o, const BaseCount& base)
	{
		cout << base.count[0];
		for (int x = 1; x < 4; x++)
			cout << ' ' << base.count[x];
		return o;
	}
};

typedef map<ContigID, Sequence> ContigMap;
typedef vector<BaseCount> BaseCounts;
typedef map<ContigID, BaseCounts> BaseCountsMap;


static ContigMap g_contigs;
static BaseCountsMap g_baseCounts;

//static const

static void readContigs(const string& contigsPath)
{
	FastaReader contigsFile(contigsPath.c_str());
	int count = 0;
	while(contigsFile.isGood()) {
		ContigID id;
		Sequence contig = contigsFile.ReadSequence(id);
		g_contigs[id] = contig;
		unsigned numBases = opt::csToNt ? contig.length() + 1 : contig.length();
		g_baseCounts[id] = BaseCounts(numBases);
		if (count == 0) {
			// Detect colour-space contigs.
			opt::colourSpace = isdigit(contig[0]);
		} else {
			if (opt::colourSpace)
				assert(isdigit(contig[0]));
			else
				assert(isalpha(contig[0]));
		}
		count++;
	}
	cerr << "size of g_contigs: " << count << '\n';
}

static int charToInt(char c)
{
	assert('0' <= c && c <= '9');
	return c - '0';
}

static void readAlignment(string& line, string& readID,
		Sequence& seq, AlignmentVector& alignments)
{
	char anchor;
	stringstream s(line);

	if (opt::colourSpace)
		s >> readID >> anchor >> seq;
	else if (opt::csToNt) {
		s >> readID >> anchor >> seq;
		seq = colourToNucleotideSpace(anchor, seq);
	} else
		s >> readID >> seq;

	Alignment alignment;
	while (s >> alignment)
		alignments.push_back(alignment);
}

static void buildBaseQuality()
{
	if (opt::csToNt)
		opt::colourSpace = false;

	for (string line; getline(cin, line);) {
		string readID;
		Sequence seq;
		AlignmentVector alignments;
		readAlignment(line, readID, seq, alignments);

		for (AlignmentVector::const_iterator alignIter = alignments.begin();
				alignIter != alignments.end(); ++alignIter) {
			const char* s;
			Alignment a;
			if (alignIter->isRC) {
				s = reverseComplement(seq).c_str();
				a = alignIter->flipQuery();
			} else {
				s = seq.c_str();
				a = *alignIter;
			}

			int read_min;
			int read_max;
			if (!opt::csToNt) {
				read_min = a.read_start_pos - a.contig_start_pos;
				read_min = read_min > 0 ? read_min : 0;

				read_max = a.read_start_pos + g_contigs[a.contig].length() -
					a.contig_start_pos;
				read_max = read_max < a.read_length ? read_max : a.read_length;
			} else {
				read_min = a.read_start_pos;
				read_max = read_min + a.align_length + 1;
				if (read_min != 0)
					continue;
			}

			for (int x = read_min; x < read_max; x++) {
				int base;
				if (!opt::colourSpace)
					base = baseToCode(s[x]);
				else
					base = charToInt(s[x]);
				unsigned loc = a.contig_start_pos - a.read_start_pos + x;
				g_baseCounts[a.contig][loc].count[base]++;
			}
		}
	}
}

static char selectBase(const BaseCount& count)
{
	int bestBase = -1;
	unsigned bestCount = 0;
	for (int x = 0; x < 4; x++)
		if (count.count[x] > bestCount) {
			bestBase = x;
			bestCount = count.count[x];
		}

	if (bestBase == -1)
		return '?';
	return codeToBase(bestBase);
}

static void consensus(const char* outPath)
{
	FastaWriter outFile(outPath);

	for (BaseCountsMap::const_iterator it = g_baseCounts.begin();
			it != g_baseCounts.end(); ++it) {
		unsigned seqLength = it->second.size();

		char outSeq[seqLength];
		memset(outSeq, '?', seqLength);

		transform(it->second.begin(), it->second.end(), outSeq, selectBase);

		LinearNumKey idKey = convertContigIDToLinearNumKey(it->first);
		outFile.WriteSequence(Sequence(outSeq, seqLength), idKey, 0.0f, it->first);
		for (unsigned i = 0; i < seqLength; i++)
			cout << idKey << ' ' << seqLength << ' ' << i
				<< ' ' << outSeq[i]
				<< ' ' << g_contigs[it->first].at(i)
				<< ' ' << g_baseCounts[it->first][i] << '\n';
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
			case 'v': opt::verbose++; break;
			case 'o': arg >> opt::outPath; break;
			case 'C': opt::csToNt = true; break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	string contigsPath(argv[argc - 1]);

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	readContigs(contigsPath);
	buildBaseQuality();

	if (opt::outPath.empty())
		consensus("consensus.fa");
	else
		consensus(opt::outPath.c_str());
}
