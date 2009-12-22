#ifndef CONTIGNODE_H
#define CONTIGNODE_H 1

#include <string>
#include <ostream>

/** A tuple of a contig ID and an orientation. */
class ContigNode {
  public:
	ContigNode() { }
	ContigNode(unsigned id, bool sense) : m_node(id << 1 | sense) { };

	unsigned id() const { return m_node >> 1; }
	bool sense() const { return m_node & 1; }

  private:
	unsigned m_node;
};

#endif
