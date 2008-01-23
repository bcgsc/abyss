#ifndef PACKEDSEQ_H
#define PACKEDSEQ_H

#include "Sequence.h"

class PackedSeq
{
	public:
		
		// Constructor/Destructor
		PackedSeq(const Sequence& seq);
		PackedSeq(char* const pData, int length);
		
		// Copy constructor
		PackedSeq(const PackedSeq& pseq);
		
		// Destructor, frees memory
		~PackedSeq();
		
		// Allocate memory for the string
		void allocate(int length);
		
		// Assignment Operator
		PackedSeq& operator=(const PackedSeq& other);
		
		// Operators
		bool operator==(const PackedSeq& other) const;
		bool operator!=(const PackedSeq& other) const;
		bool operator<(const PackedSeq& other) const;
		
		// Decode the sequence
		Sequence decode() const;
		
		// get the length of the sequence
		int getSequenceLength() const;
		
		// Get the number of bytes in the sequence
		static int getNumCodingBytes(int seqLength);
		
		// Return the pointer to the data
		const char* const getDataPtr() const;
		
		// get a particular base
		char getFirstBase() const { return getBase(0); }
		char getLastBase() const { return getBase(m_length - 1); }
		char getBase(int seqIndex) const;
		
		// Reverse and complement this sequence
		void reverseComplement();
		
		// append/prepend
		// these functions preserve the length of the sequence by shifting first before adding the new base
		// the base shifted off is returned
		char shiftAppend(char base);
		char shiftPrepend(char base);

		
		// Print
		void print() const;
		
	private:
	
		PackedSeq();
		
		// get/set a particular value
		inline void setBase(char* pSeq, int byteNum, int index, char base);
		inline char getBase(const char* pSeq, int byteNum, int index) const;
		
		// Create the two bit code for the base
		inline char baseToCode(char base) const;
		inline char codeToBase(char code) const;
		
		// complement a base
		inline int seqIndexToByteNumber(int seqIndex) const;
		inline int seqIndexToBaseIndex(int seqIndex) const;
		
		// shift a single byte
		char leftShiftByte(char* pSeq, int byteNum, int index, char base);
		char rightShiftByte(char* pSeq, int byteNum, int index, char base);
		
		// sequence is terminated by a null byte (all zeros)
		char* m_pSeq;
		char m_length;
		
};

// Global function to make a reverse complement of a packed seq
PackedSeq reverseComplement(const PackedSeq& seq);


#endif
