#ifndef CONTIGNODE_H
#define CONTIGNODE_H 1

#include "Dictionary.h"
#include "StringUtil.h"
#include <cassert>
#include <string>
#include <istream>
#include <ostream>

/** A dictionary of contig IDs. */
extern Dictionary g_contigIDs;

/** A tuple of a contig ID and an orientation. */
class ContigNode {
  public:
	ContigNode() { }
	ContigNode(unsigned id, bool sense) : m_node(id << 1 | sense) { };

	unsigned id() const { return m_node >> 1; }
	bool sense() const { return m_node & 1; }

	friend std::istream& operator >>(std::istream& in,
			ContigNode& o)
	{
		std::string s;
		if (in >> s) {
			char c = chop(s);
			assert(c == '+' || c == '-');
			o = ContigNode(g_contigIDs.serial(s), c == '-');
		}
		return in;
	}

	friend std::ostream& operator <<(std::ostream& out,
			const ContigNode& o)
	{
		return out << g_contigIDs.key(o.id())
			<< (o.sense() ? '-' : '+');
	}

  private:
	unsigned m_node;
};

#endif
