#ifndef PACKEDSEQ_H
#define PACKEDSEQ_H

#include <list>
#include "CommonUtils.h"
#include "Sequence.h"
#include "SeqExt.h"

enum SeqFlag
{
	SF_SEEN = 0x1,
	SF_DELETE = 0x2
};

typedef std::pair<int, short> SeqID;
typedef std::vector<SeqID> IDList;

class PackedSeq
{
	public:
		
		// Constructor/Destructor
		PackedSeq();
		PackedSeq(const Sequence& seq);
		PackedSeq(const char* const pData, int length);
		
		// Copy constructor
		PackedSeq(const PackedSeq& pseq);
		
		// Destructor
		~PackedSeq();
		
		// Assignment Operator
		PackedSeq& operator=(const PackedSeq& other);
		
		// Operators
		inline bool operator==(const PackedSeq& other) const
		{
			return compare(other) == 0;
		}
		
		inline bool operator!=(const PackedSeq& other) const
		{
			return compare(other) != 0;
		}
		inline bool operator<(const PackedSeq& other) const
		{
			return compare(other) < 0;	
		}
		
		// Comparison
		int compare(const PackedSeq& other) const;
		
		// Decode the sequence
		Sequence decode() const;
		Sequence decodeByte(const char byte) const;
		
		unsigned int getCode() const;
		size_t getHashCode() const;
		
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
		
		// Set a flag indicating this sequence can by extending by base b
		void setExtension(extDirection dir, SeqExt extension);
		void clearExtension(extDirection dir, char b);
		void clearAllExtensions(extDirection dir);
		bool checkExtension(extDirection dir, char b) const;
		bool hasExtension(extDirection dir) const;
		bool isAmbiguous(extDirection dir) const;
		void printExtension() const;
		
		// Reverse and complement this sequence
		void reverseComplement();
		
		// append/prepend
		// these functions preserve the length of the sequence by shifting first before adding the new base
		// the base shifted off is returned
		char rotate(extDirection dir, char base);
		char shiftAppend(char base);
		char shiftPrepend(char base);
		
		// add an ID to the list
		void addID(SeqID id);
		void addIDList(const IDList& ids);
		const IDList getIDList() const;
		
		// Print
		void print() const;
		
		// The maximum kmer size is hardcoded to be 64
		// Why is this? If we use a dynamically allocated character buffer malloc/new will give us 16 or 32 bytes no matter how much we want
		// This padding + the size of the pointer effectively negates the gains from use a compressed sequence
		// By hardcoding this value we can keep things aligned, plus remove the need for alloc/frees
		// The alternatives are a) accepting the inefficiency of small dynamic allocations or b) writing a custom small object allocator
		static const int MAX_KMER = 64;
		static const int NUM_BYTES = MAX_KMER / 4;

	private:
			
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
		
		char m_seq[NUM_BYTES];
		char m_length;
		char m_flags;
		SeqExt m_extensions[2]; // single byte each
		
#ifndef SAVE_MEM
		IDList m_ids;
#endif
};

// Global function to make a reverse complement of a packed seq
PackedSeq reverseComplement(const PackedSeq& seq);


#endif
