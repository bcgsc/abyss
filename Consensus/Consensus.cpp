#include "Aligner.h"
#include "FastaReader.h"
#include "PairUtils.h"
#include "FastaWriter.h"
#include <getopt.h>
#include <string>
#include <sstream>
#include <iostream>

using namespace std;

#define PROGRAM "Consensus"

static const char *VERSION_MESSAGE =
PROGRAM " (ABySS) " VERSION "\n"
"Written by Shaun Jackman and Tony Raymond.\n"
"\n"
"Copyright 2009 Canada's Michael Smith Genome Science Centre\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION]... [FILE]...\n"
"\n"
"Alignments and read sequences from KAligner are read in from standard\n"
"input. Ensure that the --seq option was used when running KAligner.\n"
"Write the consensus results of all reads to OUTPUT. Call a consensus\n"
"at each position of each contig and write the result to standard output.\n"
"\n"
"  -o, --out=OUTPUT      write converted sequences in fasta format to this file\n"
"  -C, --cstont          convert the colour space input to nucleotide space\n"
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
	BaseCount() { fill(count, count + 4, 0); }
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

/** Read all contigs in and store the contigs in g_contigs and make a
 * g_baseCounts, to store pile-up for each base. */
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

/** Builds the pile up of all reads based on the alignments and
 * read sequence */
static void buildBaseQuality()
{
	if (opt::csToNt)
		opt::colourSpace = false;

	// for each read and/or set of alignments.
	for (string line; getline(cin, line);) {
		string readID;
		Sequence seq;
		AlignmentVector alignments;

		readAlignment(line, readID, seq, alignments);

		// If converting to NT space, check that at least one of the
		// alignments starts at read location 0. Otherwise, it is
		// likely to introduce a frameshift or erroneous sequence in
		// the final consensus.
		if (opt::csToNt) {
			bool good = false;
			for (AlignmentVector::const_iterator alignIter = alignments.begin();
					alignIter != alignments.end(); ++alignIter) {
				if (alignIter->read_start_pos == 0) {
					good = true;
					break;
				}
			}
			if (!good)
				continue;
		}

		// For each alignment for the read.
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

			if (g_baseCounts.find(a.contig) == g_baseCounts.end())
				continue;

			BaseCounts& countsVec = g_baseCounts[a.contig];
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
			}

			if ((int)g_baseCounts[a.contig].size() < a.contig_start_pos
					- a.read_start_pos + read_max - 1)
				cerr << countsVec.size() << '\n';

			// Assertions to make sure alignment math was done right.
			assert((int)countsVec.size() >= a.contig_start_pos
					- a.read_start_pos + read_max - 1);
			assert(read_max <= (int)seq.length());
			assert(read_min >= 0);

			// Pile-up every base in the read to the contig. 
			for (int x = read_min; x < read_max; x++) {
				int base;
				if (!opt::colourSpace)
					base = baseToCode(s[x]);
				else
					base = baseToCode(s[x]);
				unsigned loc = a.contig_start_pos - a.read_start_pos + x;
				assert(loc < countsVec.size());
				countsVec[loc].count[base]++;
			}
		}
	}
}

/** Returns the most likely base found by the pile up count. */
static char selectBase(const BaseCount& count, unsigned& sumBest,
		unsigned& sumSecond)
{
	int bestBase = -1;
	unsigned bestCount = 0;
	unsigned secondCount = 0;
	for (int x = 0; x < 4; x++) {
		if (count.count[x] > bestCount) {
			bestBase = x;
			secondCount = bestCount;
			bestCount = count.count[x];
		}
	}

	sumBest += bestCount;
	sumSecond += secondCount;

	if (bestBase == -1)
		return 'N';
	return codeToBase(bestBase);
}

/** Convert all 'N' bases to nt's based on local information. */
static void fixUnknown(Sequence& ntSeq, const Sequence& csSeq )
{
	size_t index = ntSeq.find_first_of('N');
	size_t rindex = ntSeq.find_last_of('N');
	char base;
	/*if (index == 0) {
		//for (index = ntSeq.find_first_of("ACGT"); index > 0; index--)
		index = ntSeq.find_first_of("ACGT");
		while (index != 0) {
			base = colourToNucleotideSpace(ntSeq.at(index),
					csSeq.at(index - 1));
			ntSeq.replace(index - 1, 1, 1, base);
			//ntSeq[index-1] = base;
			index = ntSeq.find_first_of("ACGT");
		}
		index = ntSeq.find_first_of('N');
	}*/

	if (index == 0 || rindex == ntSeq.length() - 1) {
		ntSeq = ntSeq.substr(ntSeq.find_first_of("ACGT"),
				ntSeq.find_last_of("ACGT") -
				ntSeq.find_first_of("ACGT") + 1);
		index = ntSeq.find_first_of('N');
	}

	while (index != string::npos) {
		// If the base isn't the first or last base in the seq...
		base = colourToNucleotideSpace(ntSeq.at(index - 1),
				csSeq.at(index - 1));
		ntSeq.replace(index, 1, 1, base);
		index = ntSeq.find_first_of('N');
	}
}

/** Forms contigs based on the consensus of each base and outputs them
 * to the file specified by the -o option. */
static void consensus(const char* outPath)
{
	FastaWriter outFile(outPath);

	unsigned numIgnored = 0;
	for (BaseCountsMap::const_iterator it = g_baseCounts.begin();
			it != g_baseCounts.end(); ++it) {
		unsigned seqLength = it->second.size();

		char outSeq[seqLength];
		memset(outSeq, 'N', seqLength);

		unsigned sumBest = 0;
		unsigned sumSecond = 0;
		for (unsigned x = 0; x < seqLength; x++) {
			outSeq[x] = selectBase(it->second[x], sumBest, sumSecond);
		}

		LinearNumKey idKey = convertContigIDToLinearNumKey(it->first);
		Sequence outString = Sequence(outSeq, seqLength);

		if (outString.find_first_of("ACGT") != string::npos) {
			// Check that the average percent agreement was enough to
			// write the contig to file.
			float percentAgreement = sumBest / (float)(sumBest + sumSecond);
			if (isnan(percentAgreement) || percentAgreement < .9) {
				numIgnored++;
				if (opt::csToNt) {
					if (opt::verbose > 0)
						cerr << "Warning: Contig " << it->first
							<< " has less than 90\% agreement\n"
							<< "and will not be converted.\n";
				} else
					continue;
			} else {
				if (opt::csToNt)
					fixUnknown(outString, g_contigs[it->first]);
				outFile.WriteSequence(outString, idKey, 0);
			}

			// <contig id> <length> <consensus result> <expected> <A> <C> <G> <T>
			if (opt::csToNt)
				for (unsigned i = 0; i < seqLength - 1; i++)
					cout << idKey << ' ' << seqLength << ' ' << i << ' '
						<< nucleotideToColourSpace(outSeq[i], outSeq[i + 1])
						<< ' ' << g_contigs[it->first].at(i)
						<< ' ' << g_baseCounts[it->first][i] << '\n';
			else
				for (unsigned i = 0; i < seqLength; i++)
					cout << idKey << ' ' << seqLength << ' ' << i
						<< ' ' << outSeq[i]
						<< ' ' << g_contigs[it->first].at(i)
						<< ' ' << g_baseCounts[it->first][i] << '\n';
		} else if (opt::verbose > 0) {
			cerr << "Warning: Contig " << it->first
				<< " was not supported\n"
				<< "by a complete read and was ommited.\n";
		}
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
