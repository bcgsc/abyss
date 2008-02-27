#ifndef SEQEXT_H
#define SEQEXT_H

#include "CommonDefs.h"
#include "CommonUtils.h"

class SeqExt
{
	public:
		
		SeqExt();
		
		// Set a particular base as being present
		void SetBase(char base);
		
		// Clear a base
		void ClearBase(char base);
		
		// Check whether a base is set
		bool CheckBase(char base) const;
		
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
		unsigned char base2Bit(char base) const;
		char SeqExt::bit2Base(unsigned char code);
		
		unsigned char m_record;

	
};

#endif
