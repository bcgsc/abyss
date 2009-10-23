#ifndef Kmer_H
#define Kmer_H 1

#include "config.h"
#include "Sense.h"
#include "Sequence.h"
#include <stdint.h>

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

	unsigned getSequenceLength() const { return m_length; }

	void reverseComplement();

	bool isPalindrome() const;
	bool isPalindrome(extDirection dir) const;
	void setLastBase(extDirection dir, uint8_t base);
	uint8_t getLastBaseChar() const;

	uint8_t shift(extDirection dir, uint8_t base = 0)
	{
		return dir == SENSE ? shiftAppend(base) : shiftPrepend(base);
	}

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

#endif
