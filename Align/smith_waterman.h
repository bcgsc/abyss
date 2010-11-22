//////////////////////////////////////////////////////////////////////////////
// Smith_waterman header
//
// Written by Rong She (rshe@bcgsc.ca)
// Last modified: Jul 6, 2010
//////////////////////////////////////////////////////////////////////////////

#ifndef SMITH_WATERMAN_H
#define SMITH_WATERMAN_H

#include <string>
#include <iostream>
#include <vector>
using namespace std;

/* a simple data structure to store alignment sequences */
struct SMAlignment {
	string	query_align;
	string	target_align;
	string	match_align;

	SMAlignment() {}

	friend ostream& operator<<(ostream& os, const SMAlignment& align)
	{
		os << "query:" << align.query_align << endl;
		os << "match:" << align.match_align << endl;
		os << "targt:" << align.target_align << endl;
		return os;
	}
};

/* overlap alignment between two sequences t and h */
struct overlap_align {
	unsigned	overlap_t_pos; //overlap on t is from this pos to end of sequence t
	unsigned	overlap_h_pos; //overlap on h is from beginning of sequence h to this pos
	string		overlap_str;
	unsigned	overlap_match;

	overlap_align() : overlap_t_pos(0), overlap_h_pos(0), overlap_str(""), overlap_match(0) {}
	overlap_align(unsigned t_pos, unsigned h_pos, string& overlap, unsigned num_of_match) :
		overlap_t_pos(t_pos), overlap_h_pos(h_pos), overlap_str(overlap), overlap_match(num_of_match) {}

	unsigned length() const { return overlap_str.length(); }
	double pid() const { return (double)overlap_match / overlap_str.length(); }

	friend ostream& operator<<(ostream& os, const overlap_align& align) {
		os << "overlap region:" << align.overlap_str << endl
			<< "t:" << align.overlap_t_pos << ", h:" << align.overlap_h_pos
			<< ";num_of_match:" << align.overlap_match
			<< ";len:" << align.length() << ";pid:" << align.pid() << endl;
		return os;
	}
};

void alignOverlap(const string& seq_a, const string& seq_b,
	unsigned seq_a_start_pos, vector<overlap_align>& overlaps, bool multi_align, bool verbose);

#endif /* SMITH_WATERMAN_H */
