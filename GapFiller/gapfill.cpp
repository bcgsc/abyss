#include "dialign.h" //this has to be first...
#include "SAM.h"
#include "DataLayer/Options.h"
#include "smith_waterman.h"
#include "Align/Options.h"
#include "Uncompress.h"
#include "FastaReader.h"

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
" DIALIGN-TX options:\n"
"  -D, --dialign-d=N     dialign debug level, default: 0\n"
"  -M, --dialign-m=FILE  score matrix, default: dna_matrix.scr\n"
"  -P, --dialign-p=FILE  diagonal length probability distribution\n"
"                        default: dna_diag_prob_100_exp_550000\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

//TODO
namespace opt {
	static float identity = 0.9;
	static unsigned min_matches = 50;
	static unsigned max_overlap = 500;
	static unsigned min_size = 500;
	//static unsigned num_mismatch = 5;
	static string alignPath;

	static int dialign_debug;
	static string dialign_score;
	static string dialign_prob;

	static int verbose;
}

//TODO
static const char shortopts[] = "h:l:s:vD:M:P:";

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
	{ "dialign-d",   required_argument, NULL, 'D' },
	{ "dialign-m",   required_argument, NULL, 'M' },
	{ "dialign-p",   required_argument, NULL, 'P' },
	{ NULL, 0, NULL, 0 }
};

struct Scaffold {
	typedef pair<size_t, size_t> Gap;
	typedef vector<Gap> Gaps;

	FastaRecord rec;
	Gaps gaps;

	Scaffold(FastaRecord rec) : rec(rec) {
		splitScaffolds();
	}

	void splitScaffolds() {
		size_t j = 0;
		for (size_t i = rec.seq.find_first_of('N'); i != string::npos;
				i = rec.seq.find_first_of('N', j)) {
			j = rec.seq.find_first_not_of('N', i);
			gaps.push_back(Gap(i, j));
		}
	}
};

typedef map<string, Scaffold> Scaffolds;
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
		if (rec.seq.size() >= opt::min_size
				&& rec.seq.find_first_of("N") != string::npos) {
			stats.scaffolds++;
			assert(scaffs.find(rec.id) == scaffs.end());
			scaffs.insert(make_pair(rec.id, Scaffold(rec)));
		}
	}
}

static bool isNearGap(Scaffold::Gap& gap, SAMRecord& align)
{
	int align_start = align.pos;
	int gap_start = gap.first;
	return align_start <= gap_start && align_start >= (int)(gap_start -
				opt::max_overlap + opt::min_matches);
}

static bool isNearGaps(Scaffold& scaff, SAMRecord& align)
{
	for (Scaffold::Gaps::iterator it = scaff.gaps.begin();
			it != scaff.gaps.end(); it++)
		if (isNearGap(*it, align)) {
			return true;
		}
	return false;
}

static void readAlignments(string path, Alignments& aligns,
		Scaffolds& scaffs)
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
				&& rec.tags.find(rec.rname, 5) != string::npos
				&& scaffs.count(rec.rname) > 0) {
			stats.split_same++;
			if (isNearGaps(scaffs.find(rec.rname)->second, rec))
				aligns.insert(make_pair(rec.rname, rec));
		}
		assert(in);
	}
	assert(in.eof());
}
#if 0
static void splitScaffolds(string scaffold, vector<string>& contigs)
{
	size_t j = 0;
	for (size_t i = scaffold.find_first_of('N'); i != string::npos;
			i = scaffold.find_first_of('N', j)) {
		contigs.push_back(scaffold.substr(j, i - j));
		j = scaffold.find_first_not_of('N', i);
	}
	contigs.push_back(scaffold.substr(j));
}
bool isGapless(overlap_align& o, Sequence& s) {
	return o.length() == s.length() - o.overlap_t_pos &&
		o.length() == o.overlap_h_pos + 1;
}
#endif
static void filterGapAlignments(vector<overlap_align>& overlaps,
		string /*seq*/)
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
	if (overlaps.empty()) {
		//stats.pid_low++;
		return;
	}
#if 0
	for (it = overlaps.begin(); it != overlaps.end(); it++ ) {
		overlap_align o = *it;
		if (!isGapless(o, seq))
			overlaps.erase(it--);
	}
	if (overlaps.empty()) {
		//stats.has_indel++;
		return;
	}
#endif
}


/* @return length of overlapping read ends. <0,0> if no good alignment
*/
static void alignReadToGapFlanks(string seg1,
		string seg2, string read, vector<string>& seqs, size_t scaffold_start)
{
	//overlap align end of first segment to start of read
	vector<overlap_align> overlaps1;
	unsigned start_pos = 0;
	{
		string a(seg1), b(read);
		if (seg1.size() > opt::max_overlap)
			a = seg1.substr(seg1.size() - opt::max_overlap);
		if (read.size() > opt::max_overlap)
			b = read.substr(0, opt::max_overlap);
		alignOverlap(a, b, 0, overlaps1,
				true, opt::verbose > 2);
	}

	filterGapAlignments(overlaps1, seg1);

	//overlap align start of second segment to end of read
	vector<overlap_align> overlaps2;
	start_pos = 0;
	{
		string a(read), b(seg2);
		if (read.size() > opt::max_overlap)
			a = read.substr(read.size() - opt::max_overlap);
		if (seg2.size() > opt::max_overlap)
			b = seg2.substr(0, opt::max_overlap);
		alignOverlap(a, b, 0, overlaps2,
				true, opt::verbose > 2);
	}

	filterGapAlignments(overlaps2, seg2);

	//if both alignments have sufficient identity, return overlaps
	if (overlaps1.size() == 1 && overlaps2.size() == 1) {
		unsigned start = overlaps1[0].overlap_str.size();
		int length = read.size() - overlaps2[0].overlap_str.size() - start;
		if (length <= 0) //TODO: Handle this properly!!!
			return;
		cerr << scaffold_start << '\n';
		//cout << seg1 << '\n'
		//	<< seg2 << '\n'
		//	<< read << '\n'
		//	<< overlaps1[0] << overlaps2[0] << '\n';
		seqs.push_back(read.substr(start, length));
	}
}

static void alignReadsToGapFlanks(Scaffold& scaff,
		Alignments& aligns)
{
	string cid = scaff.rec.id;
	cerr << "examining contig " << cid << ", which has ";
	Alignments::iterator start, end;
	tie(start, end) = aligns.equal_range(cid);
	cerr << distance(start, end) << " alignments and "
		<< scaff.gaps.size() + 1 << " segments...\n";

	vector< vector<string> > gap_seqs(scaff.gaps.size());

	for (tie(start, end) = aligns.equal_range(cid); start != end;
			start++) {
		//identity the set of gaps that this alignment could span
		

		for (Scaffold::Gaps::iterator cit = scaff.gaps.begin();
				cit != scaff.gaps.end(); cit++) {
			if (!isNearGap(*cit, start->second))
				continue;
			string read_seq = start->second.seq;

			int seg1_start = max(0, (int)(cit->first - opt::max_overlap));
			assert(seg1_start >= 0);
			string seg1 = scaff.rec.seq.substr(seg1_start,
					cit->first - seg1_start);

			int seg2_end = min(scaff.rec.seq.size(),
					cit->second + opt::max_overlap);
			assert((unsigned)seg2_end <= scaff.rec.seq.size());
			string seg2 = scaff.rec.seq.substr(cit->second,
					seg2_end - cit->second);

			alignReadToGapFlanks(seg1, seg2, read_seq,
					gap_seqs[(unsigned)distance(scaff.gaps.begin(),
					cit)], cit->first);
		}
	}

	for (vector< vector<string> >::iterator it = gap_seqs.begin();
			it != gap_seqs.end(); it++) {
		string alignment;
		unsigned matches;
		if (it->size() > 0) {
			//for (vector<string>::iterator sit = it->begin();
			//		sit != it->end(); sit++) {
			//	cout << *sit << '\n';
			//	if (sit->size() <= 1)
			//		it->erase(sit--);
			//}
			//if (it->size() == 0)
			//	continue;
			Sequence consensus = dialign(*it, alignment, matches);
		//if (opt::verbose > 2)
		   	cerr << alignment << consensus << '\n';
		}
		//float identity = (float)matches / consensus.size();
		//if (identity > opt::identity) {
			//fix contig sequence
		//}
	}
}

static void fillGaps(Scaffolds& scaffs, Alignments& aligns)
{
	//foreach scaffold
	for (Scaffolds::iterator sit = scaffs.begin();
			sit != scaffs.end(); sit++) {
		//break scaffolds into contigs
		//cout << sit->second.seq << '\n';
		alignReadsToGapFlanks(sit->second, aligns);
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
			case 'D': arg >> opt::dialign_debug; break;
			case 'M': arg >> opt::dialign_score; break;
			case 'P': arg >> opt::dialign_prob; break;
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
	readAlignments(argv[optind++], aligns, scaffs);
	cerr << "aligns: " << stats.aligns << "\tsplit: "
		<< stats.split_same << '\n';

	init_parameters();
	set_parameters_dna();
	para->DEBUG = opt::dialign_debug;
	para->SCR_MATRIX_FILE_NAME = (char*)opt::dialign_score.c_str();
	para->DIAG_PROB_FILE_NAME = (char*)opt::dialign_prob.c_str();
	initDialign();

	fillGaps(scaffs, aligns);

	//fillGaps();
	//writeContigs();

	return 0;
}
