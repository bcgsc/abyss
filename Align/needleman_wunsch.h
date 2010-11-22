#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include <cassert>
#include <cctype>
#include <ostream>
#include <string>

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
		NWAlignment& align, bool verbose = false);

#endif /* NEEDLEMAN_WUNSCH_H */
