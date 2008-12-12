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

static ContigVec contigs;
static SimpleContigGraph contigGraph;

static struct {
	unsigned overlap;
	unsigned tooshort;
	unsigned homopolymer;
	unsigned motif;
	unsigned ambiguous;
} stats;

static const Sequence getSequence(LinearNumKey id, extDirection sense)
{
	const Sequence& seq = contigs[id].seq;
	return sense == SENSE ? seq : reverseComplement(seq);
}

class ContigNode {
	public:
		ContigNode(LinearNumKey id, extDirection sense)
			: m_id(id), m_sense(sense), m_seq(getSequence(id, sense))
			{}
		Sequence sequence() const { return m_seq; }
		size_t length() const { return m_seq.length(); }
		unsigned inDegree() const
		{
			return contigGraph[m_id].numEdges(!m_sense);
		}
		unsigned outDegree() const
		{
			return contigGraph[m_id].numEdges(m_sense);
		}
		string id() const
		{
			ostringstream s;
			s << m_id << (m_sense == SENSE ? '+' : '-');
			return s.str();
		}
		string idComplement() const
		{
			ostringstream s;
			s << m_id << (m_sense == SENSE ? '-' : '+');
			return s.str();
		}
	private:
		LinearNumKey m_id;
		extDirection m_sense;
		const Sequence m_seq;
};

static unsigned findOverlap(const ContigNode& t, const ContigNode& h,
		bool& mask)
{
	mask = false;
	unsigned len = min(t.length(), h.length());
	vector<unsigned> overlaps;
	overlaps.reserve(len);
	for (unsigned overlap = len; overlap >= 1; overlap--) {
		string a = t.sequence().substr(t.length()-overlap, overlap);
		string b = h.sequence().substr(0, overlap);
		if (a == b)
			overlaps.push_back(overlap);
	}
	if (overlaps.size() == 0)
		return 0;

	if (opt_verbose > 0) {
		cout << t.id() << '\t' << h.id();
		for (vector<unsigned>::const_iterator i = overlaps.begin();
				i != overlaps.end(); ++i) {
			cout << '\t' << *i;
		}
		cout << '\n';
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
#if 0
		// Allow the merge and mask the overlapping sequence.
		mask = true;
#else
		return 0;
#endif
	}

	return overlaps[0];
}

static string overlapContigs(const ContigNode& t, const ContigNode& h,
		unsigned overlap, bool mask)
{
	static unsigned n;
	unsigned id = contigs.size() + n++;

	unsigned gap = opt_k - 1 - overlap;
	string a = t.sequence().substr(t.length()-opt_k+1, gap);
	string o = h.sequence().substr(0, overlap);
	string b = h.sequence().substr(overlap, gap);
	if (mask)
		transform(o.begin(), o.end(), o.begin(), ptr_fun(::tolower));
	ostringstream s;
	s << '>' << id << ' ' << opt_k-1 + gap << " 0 "
		<< t.id() << ' ' << h.id() << " -" << overlap << '\n'
		<< a << o << b << '\n';
	return s.str();
}

static map<string, set<string> > g_edges;
static map<string, string> g_overlaps;

static void findOverlap(
		LinearNumKey refID, extDirection dir, const Estimate& est)
{
	if (est.distance >= 0)
		return;
	ContigNode ref(refID, SENSE);
	ContigNode pair(est.nID, est.isRC ? ANTISENSE : SENSE);
	const ContigNode& t = dir == SENSE ? ref : pair;
	const ContigNode& h = dir == SENSE ? pair : ref;
	unsigned overlap = 0;
	bool mask = false;
	if (t.outDegree() == 0 && h.inDegree() == 0)
		overlap = findOverlap(t, h, mask);
	if (overlap > 0) {
		g_edges[t.id()].insert(h.id());
		g_edges[h.idComplement()].insert(t.idComplement());
		if (g_overlaps.count(t.id()) == 0)
			g_overlaps[t.id()] = overlapContigs(t, h, overlap, mask);
	}
}

static string complement(const string& id)
{
	string comp(id);
	int i = comp.length() - 1;
	assert(i > 0);
	switch (comp[i]) {
		case '+': comp[i] = '-'; break;
		case '-': comp[i] = '+'; break;
		default: assert(false);
	}
	return comp;
}

static bool unambiguous(const string &t, string &h)
{
	h.clear();
	const set<string>& heads = g_edges[t];
	if (heads.size() > 1) {
		stats.ambiguous++;
		return false;
	}
	assert(heads.size() > 0);
	h = *heads.begin();
	if (g_edges[complement(h)].size() > 1) {
		stats.ambiguous++;
		return false;
	}
	return true;
}

static bool unseen(const string &t, const string &h)
{
	static map<string, string> seen;
	if (seen.count(t) > 0) {
		assert(seen[t] == h);
		return false;
	}
	string tc = complement(t);
	string hc = complement(h);
	if (seen.count(hc) > 0) {
		assert(seen[hc] == tc);
		return false;
	}
	seen.insert(make_pair(t, h));
	seen.insert(make_pair(hc, tc));
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

	while (!in.eof()) {
		EstimateRecord er;
		readEstimateRecord(in, er);
		for (size_t dirIdx = 0; dirIdx <= 1; ++dirIdx) {
			extDirection dir = (extDirection)dirIdx;
			const vector<Estimate>& ests = er.estimates[dirIdx];
			for (EstimateVector::const_iterator iter = ests.begin();
					iter != ests.end(); ++iter)
				findOverlap(er.refID, dir, *iter);
		}
	}

	for (map<string, string>::const_iterator i
			= g_overlaps.begin(); i != g_overlaps.end(); ++i) {
		const string& t = i->first;
		string h;
		if (unambiguous(t, h) && unseen(t, h)) {
			stats.overlap++;
			out << i->second;
		}
	}

	cout << "Overlap: " << stats.overlap << "\n"
		"Insignificant (<" << MINIMUM_OVERLAP << "bp): "
		<< stats.tooshort << "\n"
		"Homopolymer: " << stats.homopolymer << "\n"
		"Motif: " << stats.motif << "\n"
		"Ambiguous: " << stats.ambiguous << "\n";
	return 0;
}
