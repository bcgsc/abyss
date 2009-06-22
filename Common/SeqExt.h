#ifndef SEQEXT_H
#define SEQEXT_H 1

#include <stdint.h>

static const int NUM_BASES = 4;

static inline uint8_t complementBaseCode(uint8_t base)
{
    return ~base & 0x3;
}

class SeqExt
{
	public:
		SeqExt() : m_record(0) { };

		/** Set the specified adjacency. */
		void setBase(uint8_t base)
		{
			m_record |= 1 << base;
		}

		/** Clear the specified adjacency. */
		void clearBase(uint8_t base)
		{
			m_record &= ~(1 << base);
		}

		/** Return wheter the specified base is adjacent. */
		bool checkBase(uint8_t base) const
		{
			return m_record & (1 << base);
		}

		/** Clear all adjacency. */
		void ClearAll()
		{
			m_record = 0;
		}

		/** Return whether this kmer has any adjacent kmer. */
		bool HasExtension() const
		{
			return m_record > 0;
		}

		/** Return whether this kmer has more than one adjacent kmer.
		 */
		bool IsAmbiguous() const
		{
			bool powerOfTwo = (m_record & (m_record - 1)) > 0;
			return m_record > 0 && powerOfTwo;
		}

		void print() const;

		/** Return the complementary adjacency. */
		SeqExt complement() const;

	private:
		uint8_t m_record;
};

#endif
