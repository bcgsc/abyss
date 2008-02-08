#ifndef PACKEDSEQ_H
#define PACKEDSEQ_H

#include "Sequence.h"

enum SeqFlag
{
	SF_SEEN = 0x01,
	SF_DELETE = 0x02
};

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
		
		// get a subsequence of this packed seq
		PackedSeq subseq(int start, int len) const;		
		
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
		
		// flags
		void setFlag(SeqFlag flag);
		bool isFlagSet(SeqFlag flag) const;
		
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
		static inline void setBase(char* pSeq, int seqIndex, char base);
		static inline void setBase(char* pSeq, int byteNum, int index, char base);
		inline char getBase(const char* pSeq, int byteNum, int index) const;
		
		// Create the two bit code for the base
		static inline char baseToCode(char base);
		static inline char codeToBase(char code);
		
		// complement a base
		static inline int seqIndexToByteNumber(int seqIndex);
		static inline int seqIndexToBaseIndex(int seqIndex);
		
		// shift a single byte
		char leftShiftByte(char* pSeq, int byteNum, int index, char base);
		char rightShiftByte(char* pSeq, int byteNum, int index, char base);
		
		// sequence is terminated by a null byte (all zeros)
		char* m_pSeq;
		char m_length;
		mutable char m_flags;
		
};

// Global function to make a reverse complement of a packed seq
PackedSeq reverseComplement(const PackedSeq& seq);


#endif
