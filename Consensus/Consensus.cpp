#include "Alignment.h"
#include "Common/Options.h"
#include "ContigID.h"
#include "FastaReader.h"
#include "IOUtil.h"
#include "Uncompress.h"
#include "UnorderedMap.h"
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

using namespace std;

#define PROGRAM "Consensus"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Tony Raymond and Shaun Jackman.\n"
"\n"
"Copyright 2014 Canada's Michael Smith Genome Sciences Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... CONTIG\n"
"\n"
"Read alignments from KAligner from standard input.\n"
"Ensure that the --seq option was used when running KAligner.\n"
"Call a consensus at each position of each contig and write the\n"
"consensus in FASTA format to OUTPUT.\n"
"\n"
" Arguments:\n"
"\n"
"  CONTIG  contigs in FASTA format\n"
"\n"
" Options:\n"
"\n"
"  -o, --out=OUTPUT      write the output FASTA file to OUTPUT\n"
"  -p, --pileup=PILEUP   write the pileup to PILEUP\n"
"      --nt              output nucleotide contigs [default]\n"
"      --cs              output colour-space contigs\n"
"  -V, --variants        print only variants in the pileup\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static string outPath;
	static string pileupPath;
	static bool csToNt;
	static int outputCS;
	static int onlyVariants;
}

static const char shortopts[] = "o:p:vV";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "out",		 required_argument, NULL, 'o' },
	{ "pileup",      required_argument, NULL, 'p' },
	{ "variants",    no_argument,		&opt::onlyVariants, 1 },
	{ "nt",			 no_argument,		&opt::outputCS, 0 },
	{ "cs",			 no_argument,		&opt::outputCS, 1 },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

struct BaseCount {
	unsigned count[4];
	BaseCount() { fill(count, count + 4, 0); }

	/** Return the number of reads at this position. */
	unsigned sum() const { return accumulate(count, count+4, 0); }

	friend ostream& operator <<(ostream& out, const BaseCount& base)
	{
		out << base.count[0];
		for (int x = 1; x < 4; x++)
			out << '\t' << base.count[x];
		return out;
	}
};

typedef vector<BaseCount> BaseCounts;

struct ContigCount {
	Sequence seq;
	unsigned coverage;
	string comment;
	BaseCounts counts;
};

/** A map of contigs. The alignments reference the contig by name. */
typedef unordered_map<string, ContigCount> ContigMap;
static ContigMap g_contigs;

/** Read all contigs in and store the contigs in g_contigs and make a
 * g_baseCounts, to store pile-up for each base. */
static void readContigs(const string& contigsPath)
{
	FastaReader contigsFile(contigsPath.c_str(),
			FastaReader::NO_FOLD_CASE);
	int count = 0;
	for (FastaRecord rec; contigsFile >> rec;) {
		const Sequence& seq = rec.seq;
		ContigCount& contig = g_contigs[rec.id];
		contig.seq = seq;

		istringstream ss(rec.comment);
		unsigned length;
		contig.coverage = 0;
		ss >> length >> contig.coverage >> ws;
		getline(ss, contig.comment);

		if (count == 0) {
			// Detect colour-space contigs.
			opt::colourSpace = isdigit(seq[0]);
			if (!opt::outputCS)
				opt::csToNt = opt::colourSpace;
			else if (!opt::colourSpace) {
				cerr << "error: Cannot convert nucleotide data to "
					"colour space.\n";
				exit(EXIT_FAILURE);
			}
		} else {
			if (opt::colourSpace)
				assert(isdigit(seq[0]));
			else
				assert(isalpha(seq[0]));
		}

		contig.counts = BaseCounts(contig.seq.length()
				+ (opt::csToNt ? 1 : 0));

		count++;
	}
	cerr << "Read " << count << " contigs\n";
	assert(contigsFile.eof());
	assert(count > 0);
}

typedef vector<Alignment> AlignmentVector;

static void readAlignment(string& line, string& readID,
		Sequence& seq, AlignmentVector& alignments)
{
	char anchor;
	istringstream s(line);

	if (opt::colourSpace || opt::csToNt)
		s >> readID >> anchor >> seq;
	else
		s >> readID >> seq;

	Alignment alignment;
	while (s >> alignment)
		alignments.push_back(alignment);

	if (!alignments.empty() && opt::csToNt
			&& seq.find_first_not_of("0123") == string::npos)
		seq = colourToNucleotideSpace(anchor, seq);
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
			for (AlignmentVector::const_iterator
					alignIter = alignments.begin();
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
		for (AlignmentVector::const_iterator
				alignIter = alignments.begin();
				alignIter != alignments.end(); ++alignIter) {
			string seqrc;
			Alignment a;
			if (alignIter->isRC) {
				seqrc = reverseComplement(seq);
				a = alignIter->flipQuery();
			} else {
				seqrc = seq;
				a = *alignIter;
			}
			const char* s = seqrc.c_str();

			ContigMap::iterator contigIt = g_contigs.find(a.contig);
			if (contigIt == g_contigs.end()) {
				cerr << "error: unexpected contig ID: `" << a.contig
					<< "'\n";
				exit(EXIT_FAILURE);
			}

			BaseCounts& countsVec = contigIt->second.counts;

			int read_min;
			int read_max;
			if (!opt::csToNt) {
				read_min = a.read_start_pos - a.contig_start_pos;
				read_min = read_min > 0 ? read_min : 0;

				read_max = a.read_start_pos + countsVec.size() -
					a.contig_start_pos;
				read_max = read_max < a.read_length
					? read_max : a.read_length;
			} else {
				read_min = a.read_start_pos;
				read_max = read_min + a.align_length + 1;
			}

			if ((int)countsVec.size() < a.contig_start_pos
					- a.read_start_pos + read_max - 1)
				cerr << countsVec.size() << '\n';

			// Assertions to make sure alignment math was done right.
			assert((int)countsVec.size() >= a.contig_start_pos
					- a.read_start_pos + read_max - 1);
			assert(read_max <= (int)seq.length());
			assert(read_min >= 0);

			// Pile-up every base in the read to the contig.
			for (int x = read_min; x < read_max; x++) {
				char c = toupper(s[x]);
				switch (c) {
				  case 'A': case 'C': case 'G': case 'T':
				  case '0': case '1': case '2': case '3':
					unsigned pos
						= a.contig_start_pos - a.read_start_pos + x;
					assert(pos < countsVec.size());
					countsVec[pos].count[baseToCode(c)]++;
				}
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
#if 0
	if (index == 0) {
#if 0
		for (index = ntSeq.find_first_of("ACGT"); index > 0; index--)
#endif
		index = ntSeq.find_first_of("ACGT");
		while (index != 0) {
			base = colourToNucleotideSpace(ntSeq.at(index),
					csSeq.at(index - 1));
			ntSeq.replace(index - 1, 1, 1, base);
			//ntSeq[index-1] = base;
			index = ntSeq.find_first_of("ACGT");
		}
		index = ntSeq.find_first_of('N');
	}
#endif

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

static void writePileup(ostream& out,
		const string &id, unsigned pos,
		char refc, char genotype,
		const BaseCount& counts)
{
	char foldrefc = toupper(refc);
	if (opt::onlyVariants && foldrefc == genotype)
		return;
	out << id << '\t' // reference sequence name
		<< 1 + pos << '\t' // reference coordinate
		<< refc << '\t' // reference base
		<< genotype << '\t' // genotype
		<< "25\t" // P(genotype is wrong)
		<< "25\t" // P(genotype is the same as the reference)
		<< "25\t" // RMS mapping quality
		<< counts.sum() << '\t'; // number of reads
	switch (foldrefc) {
	  case 'A': case 'C': case 'G': case 'T':
	  case '0': case '1': case '2': case '3': {
		uint8_t ref = baseToCode(foldrefc);
		for (int i = 0; i < 4; i++)
			if (i != ref)
				out << string(counts.count[i], codeToBase(i));
		out << string(counts.count[ref], '.');
		break;
	  }
	  default:
		for (int i = 0; i < 4; i++)
				out << string(counts.count[i], codeToBase(i));
	}
	out << '\n';
	assert(out.good());
}

/** Forms contigs based on the consensus of each base and outputs them
 * to the file specified by the -o option. */
static void consensus(const string& outPath, const string& pileupPath)
{
	ofstream outFile(outPath.c_str());
	assert_good(outFile, outPath);

	ofstream pileupFile;
	ostream& pileupOut
		= pileupPath.empty() || pileupPath == "-" ? cout
		: (pileupFile.open(pileupPath.c_str()), pileupFile);
	assert_good(pileupOut, pileupPath);

	unsigned numIgnored = 0;
	for (ContigMap::const_iterator it = g_contigs.begin();
			it != g_contigs.end(); ++it) {
		const ContigCount& contig = it->second;
		unsigned seqLength = it->second.counts.size();

		Sequence outSeq(seqLength, 'N');
		unsigned sumBest = 0;
		unsigned sumSecond = 0;
		for (unsigned x = 0; x < seqLength; x++) {
			char c = selectBase(
					it->second.counts[x], sumBest, sumSecond);
			outSeq[x] = islower(contig.seq[x]) ? tolower(c) : c;
		}

		if (outSeq.find_first_of("ACGT") != string::npos) {
			// Check that the average percent agreement was enough to
			// write the contig to file.
			float percentAgreement
				= sumBest / (float)(sumBest + sumSecond);
			if (isnan(percentAgreement) || percentAgreement < .9) {
				numIgnored++;
				if (opt::csToNt) {
					if (opt::verbose > 0)
						cerr << "warning: Contig " << it->first
							<< " has less than 90% agreement "
							"and will not be converted.\n";
				} else
					continue;
			} else {
				if (opt::csToNt)
					fixUnknown(outSeq, contig.seq);
				ostringstream comment;
				comment << outSeq.length() << ' ' << contig.coverage;
				if (!contig.comment.empty())
					comment << ' ' << contig.comment;
				outFile << FastaRecord(
						it->first, comment.str(), outSeq);
				assert(outFile.good());
			}

			if (opt::verbose > 1) {
				// ID pos reference genotype A C G T
				if (opt::csToNt)
					for (unsigned i = 0; i < seqLength - 1; i++)
						cout << it->first << '\t' << 1+i
							<< '\t' << contig.seq[i]
							<< '\t' << nucleotideToColourSpace(
									outSeq[i], outSeq[i + 1])
							<< '\t' << contig.counts[i].sum()
							<< '\t' << contig.counts[i] << '\n';
				else
					for (unsigned i = 0; i < seqLength; i++)
						cout << it->first << '\t' << 1+i
							<< '\t' << contig.seq[i]
							<< '\t' << outSeq[i]
							<< '\t' << contig.counts[i].sum()
							<< '\t' << contig.counts[i] << '\n';
			}

			if (!pileupPath.empty()) {
				if (opt::csToNt)
					for (unsigned i = 0; i < seqLength-1; i++)
						writePileup(pileupOut, it->first, i,
								contig.seq[i],
								nucleotideToColourSpace(
									outSeq[i], outSeq[i+1]),
								contig.counts[i]);
				else
					for (unsigned i = 0; i < seqLength; i++)
						writePileup(pileupOut, it->first, i,
								contig.seq[i], outSeq[i],
								contig.counts[i]);
			}
		} else if (opt::verbose > 0) {
			cerr << "warning: Contig " << it->first
				<< " was not supported by a complete read "
				"and was ommited.\n";
		}
	}
}

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'v': opt::verbose++; break;
			case 'o': arg >> opt::outPath; break;
			case 'p': arg >> opt::pileupPath; break;
			case 'V': opt::onlyVariants = true; break;
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

	if (opt::outPath.empty() && opt::pileupPath.empty()) {
		cerr << PROGRAM ": " << "missing -o,--out option\n";
		die = true;
	}

	if (argc - optind < 1) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	} else if (argc - optind > 1) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}
	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	readContigs(argv[optind++]);
	buildBaseQuality();
	consensus(opt::outPath, opt::pileupPath);
}
