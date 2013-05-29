#ifndef GAPFILL_H
#define GAPFILL_H 1

#include "FastaReader.h"

namespace opt {
	static unsigned min_matches = 50;
	static unsigned max_overlap = 500;
}

struct Scaffold {
	typedef pair<size_t, size_t> Gap;
	typedef vector<Gap> Gaps;

	FastaRecord rec;
	Gaps gaps;

	Scaffold(FastaRecord rec) : rec(rec)
	{
		splitScaffold();
	}

	void splitScaffold()
	{
		size_t j = 0;
		for (size_t i = rec.seq.find_first_of('N'); i != string::npos;
				i = rec.seq.find_first_of('N', j)) {
			j = rec.seq.find_first_not_of('N', i);
			gaps.push_back(Gap(i, j));
		}
	}

	bool hasGaps() const { return gaps.size() > 0; }

	unsigned numGaps() const { return gaps.size(); }

	unsigned numSegs() const { return numGaps() + 1; }

	static bool isNearGap(const Scaffold::Gap& gap, const SAMRecord& align)
	{
		int align_start = align.pos;
		int gap_start = gap.first;
		return align_start <= gap_start && align_start >= (int)(gap_start -
				opt::max_overlap + opt::min_matches);
	}

	bool isNearGaps(const SAMRecord& align) const
	{
		for (Gaps::const_iterator it = gaps.begin();
				it != gaps.end(); it++)
			if (isNearGap(*it, align))
				return true;
		return false;
	}

	pair<unsigned, unsigned> fillGap(unsigned i, string& seq)
	{
		assert(i < gaps.size());
		Gap& g = gaps[i];
		rec.seq.replace(g.first, g.second - g.first, seq);
		return make_pair(g.second - g.first, seq.size());
	}

	friend ostream& operator<<(ostream& out, const Scaffold& scaff)
	{
		return out << scaff.rec;
	}

};
#endif
