//////////////////////////////////////////////////////////////////////////////
// Needleman-Wunsch header
//
// Written by Rong She (rshe@bcgsc.ca)
// Last modified: Jul 30, 2010
//////////////////////////////////////////////////////////////////////////////

#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include <cassert>
#include <cctype>
#include <ostream>
#include <string>

using namespace std;

/* a simple data structure to store alignment sequences */
struct NWAlignment {
	string	query_align;
	string	target_align;
	string	match_align; //consensus sequence

	NWAlignment() {}

	unsigned size() { return match_align.length(); }
	string consensus() { return match_align; }

	friend ostream& operator<<(ostream& out, const NWAlignment& o)
	{
		const string& a = o.query_align;
		const string& b = o.target_align;
		const string& c = o.match_align;
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

unsigned GetGlobalAlignment(const string& seq_a, const string& seq_b,
	NWAlignment& align, bool verbose = false);

#endif /* NEEDLEMAN_WUNSCH_H */
