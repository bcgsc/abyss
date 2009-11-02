#ifndef Kmer_H
#define Kmer_H 1

#include "config.h"
#include "Sense.h"
#include "Sequence.h"
#include <cassert>
#include <cstring> // for memcpy
#include <stdint.h>

class Kmer
{
  public:
	Kmer() { }
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

	static unsigned length() { return s_length; }

	/** Set the length of a k-mer.
	 * This value is shared by all instances.
	 */
	static void setLength(unsigned length)
	{
		assert(length <= MAX_KMER);
		assert(s_length == 0);
		assert(s_bytes == 0);
		s_length = length;
		s_bytes = (length + 3) / 4;
	}

	void reverseComplement();

	bool isPalindrome() const;
	bool isPalindrome(extDirection dir) const;
	void setLastBase(extDirection dir, uint8_t base);
	uint8_t getLastBaseChar() const;

	uint8_t shift(extDirection dir, uint8_t base = 0)
	{
		return dir == SENSE ? shiftAppend(base) : shiftPrepend(base);
	}

	/** Return the number of bytes needed. */
	static unsigned bytes() { return s_bytes; }
	static unsigned serialSize() { return NUM_BYTES; }

	size_t serialize(void* dest) const
	{
		memcpy(dest, m_seq, sizeof m_seq);
		return sizeof m_seq;
	}

	size_t unserialize(const void* src)
	{
		memcpy(m_seq, src, sizeof m_seq);
		return sizeof m_seq;
	}

  private:
	uint8_t shiftAppend(uint8_t base);
	uint8_t shiftPrepend(uint8_t base);

	uint8_t at(unsigned i) const;
	void set(unsigned i, uint8_t base);

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
	static unsigned s_length;
	static unsigned s_bytes;

	char m_seq[NUM_BYTES];
};

/** Return the reverse complement of the specified sequence. */
static inline Kmer reverseComplement(const Kmer& seq)
{
	Kmer rc(seq);
	rc.reverseComplement();
	return rc;
}

struct hashKmer
{
	size_t operator()(const Kmer& o) const { return o.getHashCode(); }
};

#endif
