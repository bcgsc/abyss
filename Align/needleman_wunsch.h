//////////////////////////////////////////////////////////////////////////////
// Needleman-Wunsch header
//
// Written by Rong She (rshe@bcgsc.ca)
// Last modified: Jul 30, 2010
//////////////////////////////////////////////////////////////////////////////

#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include <string>
#include <iostream>
#include <vector>
using namespace std;

/* a simple data structure to store alignment sequences */
struct NWAlignment {
	string	query_align;
	string	target_align;
	string	match_align; //consensus sequence

	NWAlignment() {}

	unsigned size() { return match_align.length(); }
	string consensus() { return match_align; }

	friend ostream& operator<<(ostream& os, const NWAlignment& align)
	{
		os << "query:" << align.query_align << endl;
		os << "match:" << align.match_align << endl;
		os << "targt:" << align.target_align << endl;
		return os;
	}
};

unsigned GetGlobalAlignment(const string& seq_a, const string& seq_b,
	NWAlignment& align, bool verbose = false);

#endif /* NEEDLEMAN_WUNSCH_H */
