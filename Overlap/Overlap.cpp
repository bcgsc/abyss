/**
 * Find contigs that overlap and end due to a lack of coverage.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */

#include "PairedAlgorithms.h"
#include "PairUtils.h"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

static const unsigned MINIMUM_OVERLAP = 5;

static unsigned opt_k;
static unsigned opt_verbose;
static bool opt_mask;
static bool opt_scaffold;

static ContigVec contigs;
static SimpleContigGraph contigGraph;

static struct {
	unsigned overlap;
	unsigned scaffold;
	unsigned none;
	unsigned tooshort;
	unsigned homopolymer;
	unsigned motif;
	unsigned ambiguous;
} stats;

class ContigNode
{
	public:
		ContigNode(LinearNumKey id, extDirection sense)
			: m_id(id), m_sense(sense) { }

		bool operator==(const ContigNode& b) const
		{
			return m_id == b.m_id && m_sense == b.m_sense;
		}
		bool operator<(const ContigNode b) const
		{
			return m_id != b.m_id ? m_id < b.m_id
				: m_sense < b.m_sense;
		}

		const ContigNode operator~() const
		{
			return ContigNode(m_id, !m_sense);
		}

		unsigned outDegree() const
		{
			return contigGraph[m_id].numEdges(m_sense);
		}
		unsigned inDegree() const
		{
			return contigGraph[m_id].numEdges(!m_sense);
		}

		const Sequence sequence() const
		{
			const Sequence& seq = contigs[m_id].seq;
			return m_sense == SENSE ? seq : reverseComplement(seq);
		}

		friend ostream& operator <<(ostream& o, const ContigNode a)
		{
			return o << a.m_id << (a.m_sense == SENSE ? '+' : '-');
		}

	private:
		const LinearNumKey m_id;
		const extDirection m_sense;
};

static unsigned findOverlap(const ContigNode& t_id,
		const ContigNode& h_id,
		bool& mask)
{
	mask = false;
	Sequence t = t_id.sequence();
	Sequence h = h_id.sequence();
	unsigned len = min(t.length(), h.length());
	vector<unsigned> overlaps;
	overlaps.reserve(len);
	for (unsigned overlap = len; overlap >= 1; overlap--) {
		string a = t.substr(t.length()-overlap, overlap);
		string b = h.substr(0, overlap);
		if (a == b)
			overlaps.push_back(overlap);
	}

	if (opt_verbose > 0) {
		cout << t_id << '\t' << h_id;
		for (vector<unsigned>::const_iterator i = overlaps.begin();
				i != overlaps.end(); ++i)
			cout << '\t' << *i;
		cout << '\n';
	}

	if (overlaps.size() == 0) {
		stats.none++;
		return 0;
	}

	if (overlaps[0] < MINIMUM_OVERLAP) {
		stats.tooshort++;
		return 0;
	}

	if (overlaps.size() >= 3
			&& overlaps[0]-overlaps[1] == overlaps[1]-overlaps[2]) {
		// Homopolymer run or motif.
		if (overlaps[0]-overlaps[1] == 1)
			stats.homopolymer++;
		else
			stats.motif++;
		mask = true;
	}

	return overlaps[0];
}

static string newContig(const ContigNode& t, const ContigNode& h,
		int dist, const string& seq)
{
	static unsigned n;
	unsigned id = contigs.size() + n++;
	ostringstream s;
	s << '>' << id << ' ' << seq.length() << " 0 "
		<< t << ' ' << h << ' ' << dist << '\n'
		<< seq << '\n';
	return s.str();
}

static string overlapContigs(const ContigNode& t_id,
		const ContigNode& h_id, unsigned overlap, bool mask)
{
	Sequence t = t_id.sequence();
	Sequence h = h_id.sequence();
	unsigned gap = opt_k - 1 - overlap;
	string a = t.substr(t.length()-opt_k+1, gap);
	string o = h.substr(0, overlap);
	string b = h.substr(overlap, gap);
	if (mask)
		transform(o.begin(), o.end(), o.begin(), ptr_fun(::tolower));
	return newContig(t_id, h_id, -overlap, a + o + b);
}

static string mergeContigs(const ContigNode& t, const ContigNode& h,
		const Estimate& est, unsigned overlap, bool mask)
{
	if (overlap > 0) {
		stats.overlap++;
		int dist = -overlap;
		int diff = dist - est.distance;
		if (fabs(diff) > allowedError(est.stdDev))
			cerr << "warning: overlap does not agree with estimate\n"
				<< '\t' << t << ',' << h << ' '
				<< dist << ' ' << est << endl;
		return overlapContigs(t, h, overlap, mask);
	} else if (opt_scaffold) {
		stats.scaffold++;
		if (opt_verbose > 0)
			cout << t << '\t' << h << "\t(" << est.distance << ")\n";
		string gap = est.distance <= 0 ? string("-")
			: string(est.distance, 'N');
		string ts = t.sequence();
		string hs = h.sequence();
		return newContig(t, h, est.distance,
				ts.substr(ts.length()-opt_k) + gap
				+ hs.substr(0, opt_k));
	} else
		assert(false);
}

struct Overlap {
	const ContigNode h;
	const Estimate est;
	const unsigned overlap;
	const bool mask;
	Overlap(const ContigNode& h,
			const Estimate& est, unsigned overlap, bool mask)
		: h(h), est(est), overlap(overlap), mask(mask) { }
};

static string mergeContigs(const ContigNode& t,
		const Overlap& overlap)
{
	return mergeContigs(t, overlap.h, overlap.est, overlap.overlap,
			overlap.mask);
}

static map<ContigNode, set<ContigNode> > g_edges;
static map<ContigNode, Overlap> g_overlaps;

static void findOverlap(
		LinearNumKey refID, extDirection dir, const Estimate& est)
{
	if (est.distance >= 0 && !opt_scaffold)
		return;
	ContigNode ref(refID, SENSE);
	ContigNode pair(est.nID, est.isRC ? ANTISENSE : SENSE);
	const ContigNode& t = dir == SENSE ? ref : pair;
	const ContigNode& h = dir == SENSE ? pair : ref;
	if (t.outDegree() > 0 || h.inDegree() > 0)
		return;
	bool mask;
	unsigned overlap = findOverlap(t, h, mask);
	if (mask && !opt_mask)
		return;
	if (overlap > 0 || opt_scaffold) {
		g_edges[t].insert(h);
		g_edges[~h].insert(~t);
		if (g_overlaps.count(t) == 0 && g_overlaps.count(~h) == 0)
			g_overlaps.insert(make_pair(t,
						Overlap(h, est, overlap, mask)));
	}
}

static bool unambiguous(const ContigNode &t)
{
	const set<ContigNode>& heads = g_edges[t];
	if (heads.size() > 1)
		return false;
	assert(heads.size() == 1);
	ContigNode h = *heads.begin();
	if (g_edges[~h].size() > 1)
		return false;
	assert(g_edges[~h].size() == 1);
	return true;
}

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

int main(int argc, const char *argv[])
{
	if (argc < 7) {
		cerr << "Overlap: missing arguments\n"
			"usage: Overlap K CONTIGS ADJ LEN DIST OUT\n";
		exit(EXIT_FAILURE);
	}

	if (string(argv[1]) == "-v") {
		opt_verbose++;
		argv++;
		argc--;
	}

	istringstream(argv[1]) >> opt_k;
	assert(opt_k > 0);
	string contigPath(argv[2]);
	string adjPath(argv[3]);
	string lenPath(argv[4]);
	string estPath(argv[5]);
	string outPath(argv[6]);

	PairedAlgorithms::readContigVec(contigPath, contigs);
	loadGraphFromAdjFile(&contigGraph, lenPath, adjPath);

	ofstream out(outPath.c_str());
	assert(out.is_open());
	ifstream in(estPath.c_str());
	assert_open(in, estPath);

	for (EstimateRecord er; in >> er;) {
		for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
			const vector<Estimate>& ests = er.estimates[dir];
			for (EstimateVector::const_iterator iter = ests.begin();
					iter != ests.end(); ++iter)
				findOverlap(er.refID, dir, *iter);
		}
	}
	assert(in.eof());
	in.close();

	for (map<ContigNode, Overlap>::const_iterator i
			= g_overlaps.begin(); i != g_overlaps.end(); ++i) {
		const ContigNode& t = i->first;
		if (unambiguous(t)) {
			out << mergeContigs(t, i->second);
			assert(out.good());
		} else
			stats.ambiguous++;
	}
	out.close();

	cout << "Overlap: " << stats.overlap << "\n"
		"Scaffold: " << stats.scaffold << "\n"
		"No overlap: " << stats.none << "\n"
		"Insignificant (<" << MINIMUM_OVERLAP << "bp): "
		<< stats.tooshort << "\n"
		"Homopolymer: " << stats.homopolymer << "\n"
		"Motif: " << stats.motif << "\n"
		"Ambiguous: " << stats.ambiguous << "\n";
	return 0;
}
