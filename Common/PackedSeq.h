#ifndef PACKEDSEQ_H
#define PACKEDSEQ_H

#include "config.h"
#include "Sense.h"
#include "SeqExt.h"
#include "Sequence.h"
#include <stdint.h>
#include <vector>

enum SeqFlag
{
	SF_MARK_SENSE = 0x1,
	SF_MARK_ANTISENSE = 0x2,
	SF_DELETE = 0x4,
};


//
// Sense extension is index 0, antisense extension is index 1
//
struct ExtensionRecord
{
	SeqExt dir[2];	
};

class PackedSeq
{
	public:
		PackedSeq();
		explicit PackedSeq(const Sequence& seq);

		// Write this packed sequence to the buffer
		size_t serialize(char* buffer) const;
		
		// Read this packed sequence from the buffer
		size_t unserialize(const char* buffer);		
		
		/** The size of the serialized structure. */
		static size_t serialSize() {
			PackedSeq *p = NULL;
			return sizeof p->m_seq + sizeof p->m_length
				+ sizeof p->m_flags + sizeof p->m_multiplicity
				+ sizeof p->m_extRecord;
		}

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

		unsigned getCode() const;
		size_t getHashCode() const;
		
		// get a subsequence of this packed seq
		PackedSeq subseq(unsigned start, unsigned len) const;		
		
		// get the length of the sequence
		unsigned getSequenceLength() const;
		
		// get the multiplicity of this sequence
		unsigned getMultiplicity(extDirection dir) const
		{
			return m_multiplicity[dir];
		}

		// get the multiplicity of this sequence
		unsigned getMultiplicity() const
		{
			return m_multiplicity[SENSE] + m_multiplicity[ANTISENSE];
		}

		// add to the multiplicity
		void addMultiplicity(extDirection dir)
		{
			if (m_multiplicity[dir] < 65535)
				++m_multiplicity[dir];
			assert(m_multiplicity[dir] > 0);
		}

		/** Set the multiplicity (not strand specific). */
		void setMultiplicity(unsigned multiplicity)
		{
			m_multiplicity[SENSE] = multiplicity;
			m_multiplicity[ANTISENSE] = 0;
		}

		uint8_t getBaseCode(unsigned seqIndex) const;
		uint8_t getLastBaseChar() const;

		// flags
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

		// Extension management
		SeqExt getExtension(extDirection dir) const;
		ExtensionRecord extension() const { return m_extRecord; }
		void setBaseExtension(extDirection dir, uint8_t base);
		void clearExtension(extDirection dir, uint8_t base);
		void clearAllExtensions(extDirection dir);
		bool hasExtension(extDirection dir) const;
		bool isAmbiguous(extDirection dir) const;

		// Reverse and complement this sequence
		void reverseComplement();
		
		// append/prepend
		// these functions preserve the length of the sequence by shifting first before adding the new base
		// the base shifted off is returned
		uint8_t shift(extDirection dir, uint8_t base = 0);
		uint8_t shiftAppend(uint8_t base);
		uint8_t shiftPrepend(uint8_t base);
		void setLastBase(extDirection dir, uint8_t base);

		bool isPalindrome() const;
		bool isPalindrome(extDirection dir) const;

#if MAX_KMER > 96
# error MAX_KMER must be no larger than 96.
#endif
#if MAX_KMER % 4 != 0
# error MAX_KMER must be a multiple of 4.
#endif
		static const unsigned NUM_BYTES = MAX_KMER / 4;

	private:
		// get/set a particular value
		static inline void setBaseCode(char* pSeq,
				unsigned seqIndex, uint8_t code);
		static inline void setBaseCode(char* pSeq,
				unsigned byteNum, unsigned index, uint8_t code);
		static inline uint8_t getBaseCode(const char* pSeq,
				unsigned byteNum, unsigned index);

		// Get the number of bytes in the sequence
		static inline unsigned getNumCodingBytes(unsigned seqLength);
		static inline unsigned seqIndexToByteNumber(unsigned seqIndex);
		static inline unsigned seqIndexToBaseIndex(unsigned seqIndex);

		// shift a single byte
		static uint8_t leftShiftByte(char* pSeq,
				unsigned byteNum, unsigned index, uint8_t base);
		static uint8_t rightShiftByte(char* pSeq,
				unsigned byteNum, unsigned index, uint8_t base);

		char m_seq[NUM_BYTES];
		uint8_t m_length;
		char m_flags;
		uint16_t m_multiplicity[2];
		ExtensionRecord m_extRecord;
};

// Global function to make a reverse complement of a packed seq
PackedSeq reverseComplement(const PackedSeq& seq);

typedef std::vector<PackedSeq> PSequenceVector;
typedef PSequenceVector::iterator PSequenceVectorIterator;

// Hash/Set functions
struct PackedSeqEqual
{
	bool operator()(const PackedSeq& obj1, const PackedSeq& obj2) const;	
};

struct PackedSeqHasher
{
	size_t operator()(const PackedSeq& myObj) const;
};

#endif
