#ifndef ASSEMBLY_VERTEXDATA_H
#define ASSEMBLY_VERTEXDATA_H 1

#include "Common/Sense.h"
#include <cassert>
#include <stdint.h>
#include <ostream>

enum SeqFlag
{
	SF_MARK_SENSE = 0x1,
	SF_MARK_ANTISENSE = 0x2,
	SF_DELETE = 0x4,
};

static inline SeqFlag complement(SeqFlag flag)
{
	unsigned out = 0;
	if (flag & SF_MARK_SENSE)
		out |= SF_MARK_ANTISENSE;
	if (flag & SF_MARK_ANTISENSE)
		out |= SF_MARK_SENSE;
	if (flag & SF_DELETE)
		out |= SF_DELETE;
	return SeqFlag(out);
}

/** Vertex properties of a de Bruijn Graph vertex. */
template <typename S, typename Set>
class VertexData
{
/** Maximum value of k-mer coverage. */
#define COVERAGE_MAX 32767U

  public:
	/** A symbol, such as a nucleotide or amino acid. */
	typedef S Symbol;

	/** A set of symbols. */
	typedef Set SymbolSet;

	/** A pair of edge sets; one for out edges and one for in edges. */
	struct SymbolSetPair
	{
		SymbolSet dir[2];
		SymbolSetPair complement() const
		{
			SymbolSetPair o;
			o.dir[SENSE] = dir[ANTISENSE].complement();
			o.dir[ANTISENSE] = dir[SENSE].complement();
			return o;
		}
	};

	VertexData() : m_flags(0)
	{
		m_multiplicity[SENSE] = 1;
		m_multiplicity[ANTISENSE] = 0;
	}

	VertexData(extDirection dir, unsigned multiplicity) : m_flags(0)
	{
		assert(multiplicity <= COVERAGE_MAX);
		m_multiplicity[dir] = multiplicity;
		m_multiplicity[!dir] = 0;
	}

	VertexData(unsigned multiplicity, SymbolSetPair ext)
		: m_flags(0), m_ext(ext)
	{
		setMultiplicity(multiplicity);
	}

	unsigned getMultiplicity(extDirection dir) const
	{
		return m_multiplicity[dir];
	}

	unsigned getMultiplicity() const
	{
		return m_multiplicity[SENSE] + m_multiplicity[ANTISENSE];
	}

	void addMultiplicity(extDirection dir, unsigned n = 1)
	{
		m_multiplicity[dir]
			= std::min(m_multiplicity[dir] + n, COVERAGE_MAX);
		assert(m_multiplicity[dir] > 0);
	}

	/** Set the multiplicity (not strand specific). */
	void setMultiplicity(unsigned multiplicity)
	{
		assert(multiplicity <= 2*COVERAGE_MAX);
		// Split the multiplicity over both senses.
		m_multiplicity[SENSE] = (multiplicity + 1) / 2;
		m_multiplicity[ANTISENSE] = multiplicity / 2;
		assert(getMultiplicity() == multiplicity);
	}

	void setFlag(SeqFlag flag) { m_flags |= flag; }
	bool isFlagSet(SeqFlag flag) const { return m_flags & flag; }
	void clearFlag(SeqFlag flag) { m_flags &= ~flag; }

	/** Return true if the specified sequence is deleted. */
	bool deleted() const { return isFlagSet(SF_DELETE); }

	/** Return true if the specified sequence is marked. */
	bool marked(extDirection sense) const
	{
		return isFlagSet(sense == SENSE
				? SF_MARK_SENSE : SF_MARK_ANTISENSE);
	}

	/** Return true if the specified sequence is marked. */
	bool marked() const
	{
		return isFlagSet(SeqFlag(SF_MARK_SENSE | SF_MARK_ANTISENSE));
	}

	SymbolSetPair extension() const { return m_ext; }

	SymbolSet getExtension(extDirection dir) const
	{
		return m_ext.dir[dir];
	}

	void setBaseExtension(extDirection dir, Symbol x)
	{
		m_ext.dir[dir].setBase(x);
	}

	void removeExtension(extDirection dir, SymbolSet ext)
	{
		m_ext.dir[dir].clear(ext);
	}

	bool hasExtension(extDirection dir) const
	{
		return m_ext.dir[dir].hasExtension();
	}

	bool isAmbiguous(extDirection dir) const
	{
		return m_ext.dir[dir].isAmbiguous();
	}

	/** Return the complement of this data. */
	VertexData operator~() const
	{
		VertexData o;
		o.m_flags = complement(SeqFlag(m_flags));
		o.m_multiplicity[0] = m_multiplicity[1];
		o.m_multiplicity[1] = m_multiplicity[0];
		o.m_ext = m_ext.complement();
		return o;
	}

	friend std::ostream& operator<<(
			std::ostream& out, const VertexData& o)
	{
		return out << "C="
			<< o.m_multiplicity[0] + o.m_multiplicity[1];
	}

  private:
	uint8_t m_flags;
	uint16_t m_multiplicity[2];
	SymbolSetPair m_ext;
};

#endif
