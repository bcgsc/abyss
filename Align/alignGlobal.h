#ifndef ALIGNGLOBAL_H
#define ALIGNGLOBAL_H

#include <cassert>
#include <cctype>
#include <utility>
#include <iostream>
#include <string>
#include <vector>

/** The result of a Needleman-Wunsch alignment. */
struct NWAlignment {
	std::string query_align;
	std::string target_align;
	std::string match_align; //consensus sequence

	NWAlignment() {}

	unsigned size() { return match_align.length(); }
	std::string consensus() { return match_align; }

	friend std::ostream& operator<<(std::ostream& out,
			const NWAlignment& o)
	{
		const std::string& a = o.query_align;
		const std::string& b = o.target_align;
		const std::string& c = o.match_align;
		assert(a.size() == c.size());
		assert(b.size() == c.size());
		for (unsigned i = 0; i < c.size(); ++i)
			out << (toupper(a[i]) == toupper(c[i]) ? '.' : a[i]);
		out << '\n';
		for (unsigned i = 0; i < c.size(); ++i)
			out << (toupper(b[i]) == toupper(c[i]) ? '.' : b[i]);
		out << '\n';
		return out << c << '\n';
	}
};

unsigned alignGlobal(
		const std::string& a, const std::string& b,
		NWAlignment& align);

/** Align the specified pair of sequences.
 * @return the number of matches and size of the consensus
 */
static inline std::pair<unsigned, unsigned> alignPair(
		const std::string& seqa, const std::string& seqb, NWAlignment& align)
{
	unsigned matches = alignGlobal(seqa, seqb, align);
	return std::make_pair(matches, align.size());
}

/** Align the specified sequences.
 * @return the number of matches and size of the consensus
 */
template <typename Seq>
static std::pair<unsigned, unsigned> alignMulti(
		const std::vector<Seq>& seqs, NWAlignment& align)
{
	Seq alignment = seqs[0];
	unsigned matches = 0;
	for (unsigned j = 0; j < seqs.size() - 1; j++) {
		matches = std::min(matches, alignGlobal(alignment,
					seqs[j+1], align));
		alignment = align.match_align;
	}
	return std::make_pair(matches, alignment.size());
}

/** Align the specified sequences.
 * @return the number of matches and size of the consensus
 */
template <typename Seq>
static std::pair<unsigned, unsigned> align(
		const std::vector<Seq>& seqs, NWAlignment& aln)
{
	assert(seqs.size() > 1);
	if (seqs.size() == 2)
		return alignPair(seqs[0], seqs[1], aln);
	else
		return alignMulti(seqs, aln);
}

template <typename Seq>
static std::pair<unsigned, unsigned> align(const std::vector<Seq>& seqs)
{
	NWAlignment aln;
	return align(seqs, aln);
}

#endif
