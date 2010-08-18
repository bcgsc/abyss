#ifndef KMERDATA_H
#define KMERDATA_H 1

#include "Sense.h"
#include "SeqExt.h"
#include <cassert>
#include <stdint.h>

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

struct ExtensionRecord
{
	SeqExt dir[2];
	ExtensionRecord operator ~()
	{
		ExtensionRecord o;
		o.dir[SENSE] = dir[ANTISENSE].complement();
		o.dir[ANTISENSE] = dir[SENSE].complement();
		return o;
	}
};

class KmerData
{
  public:
	KmerData() : m_flags(0)
	{
		m_multiplicity[SENSE] = 1;
		m_multiplicity[ANTISENSE] = 0;
	}

	KmerData(unsigned multiplicity, ExtensionRecord ext)
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

	static const unsigned COVERAGE_MAX = 32767;

	void addMultiplicity(extDirection dir)
	{
		if (m_multiplicity[dir] < COVERAGE_MAX)
			++m_multiplicity[dir];
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

	ExtensionRecord extension() const { return m_ext; }

	SeqExt getExtension(extDirection dir) const
	{
		return m_ext.dir[dir];
	}

	void setBaseExtension(extDirection dir, uint8_t base)
	{
		m_ext.dir[dir].setBase(base);
	}

	void removeExtension(extDirection dir, SeqExt ext)
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

  protected:
	uint8_t m_flags;
	uint16_t m_multiplicity[2];
	ExtensionRecord m_ext;
};

#endif
