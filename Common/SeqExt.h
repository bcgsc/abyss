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
		SeqExt();

		// Set a particular base as being present
		void setBase(uint8_t base);

		// Clear a base
		void clearBase(uint8_t base);

		// Check whether a base is set
		bool checkBase(uint8_t base) const;

		// Clear all the bits
		void ClearAll();
		
		// Check whether the sequence has any extension
		bool HasExtension() const;
		
		// Check whether the sequence has more than 1 extension
		bool IsAmbiguous() const; 
		
		void print() const;
		
		// Return a seqext object which is complementary to this one (swaps A/T, G/C)
		SeqExt complement() const;
		
	private:
		uint8_t m_record;
};

#endif
