#ifndef PACKEDSEQ_H
#define PACKEDSEQ_H 1

#include "config.h"
#include "Sense.h"
#include "SeqExt.h"
#include "Sequence.h"
#include <cassert>
#include <stdint.h>
#include <vector>

enum SeqFlag
{
	SF_MARK_SENSE = 0x1,
	SF_MARK_ANTISENSE = 0x2,
	SF_DELETE = 0x4,
};

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

class Kmer
{
  public:
	Kmer() : m_length(0) { }
	explicit Kmer(const Sequence& seq);

	int compare(const Kmer& other) const;

	bool operator==(const Kmer& other) const
	{
		return compare(other) == 0;
	}

	bool operator!=(const Kmer& other) const
	{
		return compare(other) != 0;
	}

	bool operator<(const Kmer& other) const
	{
		return compare(other) < 0;
	}

	Sequence decode() const;

	unsigned getCode() const;
	size_t getHashCode() const;

	Kmer subseq(unsigned start, unsigned len) const;

	unsigned getSequenceLength() const;

	void reverseComplement();

	bool isPalindrome() const;
	bool isPalindrome(extDirection dir) const;
	uint8_t shift(extDirection dir, uint8_t base = 0);
	void setLastBase(extDirection dir, uint8_t base);
	uint8_t getLastBaseChar() const;

  private:
	uint8_t shiftAppend(uint8_t base);
	uint8_t shiftPrepend(uint8_t base);

	static inline void setBaseCode(char* pSeq,
			unsigned seqIndex, uint8_t code);
	static inline void setBaseCode(char* pSeq,
			unsigned byteNum, unsigned index, uint8_t code);
	static inline uint8_t getBaseCode(const char* pSeq,
			unsigned byteNum, unsigned index);
	uint8_t getBaseCode(unsigned seqIndex) const;

	static inline unsigned getNumCodingBytes(unsigned seqLength);
	static inline unsigned seqIndexToByteNumber(unsigned seqIndex);
	static inline unsigned seqIndexToBaseIndex(unsigned seqIndex);

	static uint8_t leftShiftByte(char* pSeq,
			unsigned byteNum, unsigned index, uint8_t base);
	static uint8_t rightShiftByte(char* pSeq,
			unsigned byteNum, unsigned index, uint8_t base);

  public:
#if MAX_KMER > 96
# error MAX_KMER must be no larger than 96.
#endif
#if MAX_KMER % 4 != 0
# error MAX_KMER must be a multiple of 4.
#endif
	static const unsigned NUM_BYTES = MAX_KMER / 4;

  protected:
	char m_seq[NUM_BYTES];
	uint8_t m_length;
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
		: m_flags(0), m_extRecord(ext)
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

	void addMultiplicity(extDirection dir)
	{
		if (m_multiplicity[dir] < 65535)
			++m_multiplicity[dir];
		assert(m_multiplicity[dir] > 0);
	}

	/** Set the multiplicity (not strand specific). */
	void setMultiplicity(unsigned multiplicity)
	{
		assert(multiplicity <= 2*65535);
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
	bool marked(extDirection sense = SENSE) const
	{
		return isFlagSet(sense == SENSE
				? SF_MARK_SENSE : SF_MARK_ANTISENSE);
	}

	SeqExt getExtension(extDirection dir) const;
	ExtensionRecord extension() const { return m_extRecord; }
	void setBaseExtension(extDirection dir, uint8_t base);
	void removeExtension(extDirection dir, SeqExt ext);
	bool hasExtension(extDirection dir) const;
	bool isAmbiguous(extDirection dir) const;

  protected:
	char m_flags;
	uint16_t m_multiplicity[2];
	ExtensionRecord m_extRecord;
};

class PackedSeq : public Kmer, public KmerData
{
  public:
	PackedSeq() { }
	PackedSeq::PackedSeq(const Sequence& seq)
		: Kmer(seq) { }

	/** Create a PackedSeq of the specified key and data. */
	PackedSeq(const Kmer& key,
			unsigned multiplicity, ExtensionRecord ext)
		: Kmer(key), KmerData(multiplicity, ext)
	{
	}

	size_t serialize(char* buffer) const;
	size_t unserialize(const char* buffer);

	/** The size of the serialized structure. */
	static size_t serialSize() {
		PackedSeq *p = NULL;
		return sizeof p->m_seq + sizeof p->m_length
			+ sizeof p->m_flags + sizeof p->m_multiplicity
			+ sizeof p->m_extRecord;
	}
};

PackedSeq reverseComplement(const PackedSeq& seq);

typedef std::vector<PackedSeq> PSequenceVector;
typedef PSequenceVector::iterator PSequenceVectorIterator;

struct PackedSeqEqual
{
	bool operator()(const PackedSeq& a, const PackedSeq& b) const
	{
		return a == b;
	}
};

struct PackedSeqHasher
{
	size_t operator()(const PackedSeq& o) const
	{
		return o.getHashCode();
	}
};

#endif
