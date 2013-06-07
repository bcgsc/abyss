#include "alignGlobal.h"
#include "SAM.h"
#include "DataLayer/Options.h"
#include "smith_waterman.h"
#include "Align/Options.h"
#include "Uncompress.h"
#include "FastaReader.h"
#include "gapfill.h"

#include "config.h"
#include <cstdlib>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <boost/tuple/tuple.hpp>


using namespace std;
using namespace boost;

#define PROGRAM "abyss-gapfill"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Anthony Raymond.\n"
"\n"
"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... CONTIGS ALIGNS\n"
"Attempts to fill gaps in CONTIGS with spanning sequences\n"
"from ALIGNS.\n"
"\n"
"  -l, --min-align=N     the minimal alignment size [1]\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static float identity = 0.9;
	static unsigned min_size = 500;

	static int verbose;
}

static const char shortopts[] = "l:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "min-align", required_argument, NULL, 'l' },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};


typedef map<string, Scaffold> Scaffolds;
typedef multimap<string, SAMRecord> Alignments;

static struct {
	int seqs;
	int scaffolds;
	int aligns;
	int split_same;
	int gaps;
	int gaps_filled;
	int n_removed;
	int bases_added;
} stats;

static void readScaffolds(const char* path, Scaffolds& scaffs)
{
	FastaReader in(path, FastaReader::NO_FOLD_CASE);
	FastaRecord rec;
	if (opt::verbose)
		cerr << "Loading scaffolds from `" << path << "'...\n";

	while (in >> rec) {
		if (++stats.seqs % 10000000 == 0 && opt::verbose)
			cerr << "Loaded " << stats.scaffolds << " scaffolds "
				<< "out of " << stats.seqs << " sequences.\n";

		// Store only long scaffolded sequences
		if (rec.seq.size() >= opt::min_size) {
			assert(scaffs.find(rec.id) == scaffs.end());
			Scaffold scaff(rec);
			if (scaff.hasGaps()) {
				stats.scaffolds++;
				stats.gaps += scaff.numGaps();
				scaffs.insert(make_pair(rec.id, scaff));
			}
		}
	}
}

static void readAlignments(const char* path, Alignments& aligns,
		Scaffolds& scaffs)
{
	ifstream in(path);
	SAMRecord rec;
	if (opt::verbose)
		cerr << "Loading alignments from `" << path << "'...\n";

	// Parse out the headers.
	while (in.peek() == '@') {
		string tmp;
		getline(in, tmp);
		assert(in);
	}

	while (in >> rec) {
		if (++stats.aligns % 10000000 == 0 && opt::verbose)
			cerr << "Loaded " << stats.split_same << " alignments "
				<< "out of " << stats.aligns << ".\n";

		// If aligns to either side of gap store in aligns.
		if (rec.tags.size() >= 2
				&& rec.tags[0] == 'X'
				&& rec.tags[1] == 'A'
				&& rec.tags.find(rec.rname, 5) != string::npos
				&& scaffs.count(rec.rname) > 0
				&& scaffs.find(rec.rname)->second.isNearGaps(rec)) {
			stats.split_same++;
			aligns.insert(make_pair(rec.rname, rec));
		}
		assert(in);
	}
	assert(in.eof());
}

static void filterGapAlignments(vector<overlap_align>& overlaps)
{
	if (overlaps.empty()) {
		//stats.no_alignment++;
		return;
	}

	vector<overlap_align>::iterator it;
	for (it = overlaps.begin(); it != overlaps.end(); it++ ) {
		overlap_align o = *it;
		if (o.overlap_match < opt::min_matches)
			overlaps.erase(it--);
	}
	if (overlaps.empty()) {
		//stats.low_matches++;
		return;
	}

	for (it = overlaps.begin(); it != overlaps.end(); it++ ) {
		overlap_align o = *it;
		if (o.pid() < opt::identity)
			overlaps.erase(it--);
	}
}

static void alignReadToGapFlanks(string seg1,
		string seg2, string read, vector<string>& seqs)
{
	//overlap align end of first segment to start of read
	vector<overlap_align> overlaps1;
	{
		string a(seg1), b(read);
		if (seg1.size() > opt::max_overlap)
			a = seg1.substr(seg1.size() - opt::max_overlap);
		if (read.size() > opt::max_overlap)
			b = read.substr(0, opt::max_overlap);
		alignOverlap(a, b, 0, overlaps1,
				true, opt::verbose > 2);
	}

	filterGapAlignments(overlaps1);

	//overlap align start of second segment to end of read
	vector<overlap_align> overlaps2;
	{
		string a(read), b(seg2);
		if (read.size() > opt::max_overlap)
			a = read.substr(read.size() - opt::max_overlap);
		if (seg2.size() > opt::max_overlap)
			b = seg2.substr(0, opt::max_overlap);
		alignOverlap(a, b, 0, overlaps2,
				true, opt::verbose > 2);
	}

	filterGapAlignments(overlaps2);

	//if both alignments have sufficient identity, return overlaps
	if (overlaps1.size() == 1 && overlaps2.size() == 1) {
		unsigned start = overlaps1[0].overlap_str.size();
		int length = read.size() - overlaps2[0].overlap_str.size() - start;
		if (length <= 0) //TODO: Handle overlapping scaffolds properly!!!
			return;
		seqs.push_back(read.substr(start, length));
	}
}

static void alignReadsToGapFlanks(Scaffold& scaff,
		const Alignments& aligns)
{
	string cid = scaff.rec.id;
	Alignments::const_iterator start, end;
	tie(start, end) = aligns.equal_range(cid);
	if (opt::verbose > 1)
		cerr << "examining contig " << cid << ", which has "
			<< distance(start, end) << " alignments and "
			<< scaff.numSegs() << " segments...\n";

	vector< vector<string> > gap_seqs(scaff.gaps.size());

	//identify the set of gaps that this alignment could span
	for (unsigned i = 0; i < scaff.gaps.size(); i++) {
		vector<string>& seqs = gap_seqs[i];
		for ( ; start != end; start++) {
			Scaffold::Gap& gap = scaff.gaps[i];
			if (!Scaffold::isNearGap(gap, start->second))
				continue;
			string read_seq = start->second.seq;

			int seg1_start = max(0, (int)(gap.first - opt::max_overlap));
			assert(seg1_start >= 0);
			string seg1 = scaff.rec.seq.substr(seg1_start,
					gap.first - seg1_start);

			int seg2_end = min(scaff.rec.seq.size(),
					gap.second + opt::max_overlap);
			assert((unsigned)seg2_end <= scaff.rec.seq.size());
			string seg2 = scaff.rec.seq.substr(gap.second,
					seg2_end - gap.second);

			alignReadToGapFlanks(seg1, seg2, read_seq,
					seqs);
		}
	}

	for (int i = gap_seqs.size() - 1; i >= 0; i--) {
		vector<string>& seqs = gap_seqs[i];
		if (seqs.size() == 0)
			continue;
		string alignment;
		unsigned matches;
		switch (seqs.size()) {
			case 1:
				alignment = seqs[0];
				break;
			default:
				NWAlignment align;
				alignment = seqs[0];
				for (unsigned j = 0; j < seqs.size() - 1; j++) {
					matches = max(matches, alignGlobal(alignment,
								seqs[j+1], align));
					alignment = align.match_align;
				}
				break;
				//if using dialign:
				//  Sequence consensus = dialign(*it, alignment, matches);
		}
		if (alignment != "") {
			pair<unsigned, unsigned> result =
				scaff.fillGap(i,alignment);
			stats.n_removed += result.first;
			stats.bases_added += result.second;
			stats.gaps_filled++;
		}
	}
}

static void fillGaps(Scaffolds& scaffs, const Alignments& aligns)
{
	//foreach scaffold
	for (Scaffolds::iterator sit = scaffs.begin();
			sit != scaffs.end(); sit++) {
		alignReadsToGapFlanks(sit->second, aligns);
	}
}

static void printFixedContigs(const Scaffolds& scaffs,
		const char* path, ostream& out)
{
	FastaReader in(path, FastaReader::NO_FOLD_CASE);
	FastaRecord rec;
	cerr << "Fixing scaffolds from `" << path << "'...\n";

	while (in >> rec) {
		Scaffolds::const_iterator it = scaffs.find(rec.id);
		if (it == scaffs.end())
			out << rec;
		else
			out << it->second;
	}
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
	readScaffolds(argv[argc-2], scaffs);
	if (opt::verbose)
		cerr << "Loaded " << stats.scaffolds << " scaffolds "
			<< "out of " << stats.seqs << " sequences.\n";

	Alignments aligns;
	readAlignments(argv[argc-1], aligns, scaffs);
	if (opt::verbose)
		cerr << "aligns: " << stats.aligns << "\tsplit: "
			<< stats.split_same << '\n';

	fillGaps(scaffs, aligns);

	if (opt::verbose)
		cerr << "Contigs: " << stats.seqs
			<< "\nScaffolds: " << stats.scaffolds
			<< "\nGaps: " << stats.gaps
			<< "\nGaps filled: " << stats.gaps_filled
			<< "\nN's removed: " << stats.n_removed
			<< "\nBases added: " << stats.bases_added << '\n';

	printFixedContigs(scaffs, argv[argc-2], cout);

	return 0;
}
