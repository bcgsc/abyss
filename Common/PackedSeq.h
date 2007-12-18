#ifndef PACKEDSEQ_H
#define PACKEDSEQ_H

#include "Sequence.h"

class PackedSeq
{
	public:
		
		// Constructor/Destructor
		PackedSeq(const Sequence& seq);
		PackedSeq(char* const pData, int length);
		~PackedSeq();
		
		// Decode the sequence
		Sequence decode() const;
		
		// get the length of the sequence
		int getSequenceLength() const;
		
		// Get the number of bytes in the sequence
		static int getNumCodingBytes(int seqLength);
		
		// Return the pointer to the data
		const char* const getDataPtr() const;
		
		// Reverse and complement this sequence
		void reverseComplement();
		
		// Extend the sequence left/right by the given base
		// Returns a new sequence
		
		// Print
		void print() const;
		
	private:
		
		// get/set a particular value
		inline void setBase(char* pSeq, int byteNum, int index, char base);
		inline char getBase(const char* pSeq, int byteNum, int index) const;
		
		// Create the two bit code for the base
		inline char baseToCode(char base) const;
		inline char codeToBase(char code) const;
		
		// complement a base
		inline int seqIndexToByteNumber(int seqIndex) const;
		inline int seqIndexToBaseIndex(int seqIndex) const;
		
		


		// sequence is terminated by a null byte (all zeros)
		char* m_pSeq;
		char m_length;
		
};

#endif
