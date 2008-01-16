#include <math.h>
#include "PackedSeq.h"

// Default constructor
PackedSeq::PackedSeq() : m_pSeq(0), m_length(0)
{
	
}

// Destructor
PackedSeq::~PackedSeq()
{
	delete [] m_pSeq;
	m_pSeq = NULL;	
}

// Construct a sequence from a series of bytes
PackedSeq::PackedSeq(char* const pData, int length)
{
	// Allocate space for the sequence
	int numBytes = getNumCodingBytes(length);
	m_pSeq = new char[numBytes];
	
	// copy over the bytes
	memcpy(m_pSeq, pData, numBytes);
	
	// set the sequence length
	m_length = length;
}

// Construct a sequence from a String-based sequence
PackedSeq::PackedSeq(const Sequence& seq)
{
	int length = seq.length();
	
	assert(length < 255);
	
	//printf("packing %s\n", seq.c_str());
	// calculate the number of triplets required

	m_length = length;
	int numBytes = getNumCodingBytes(length);

	
	//printf("length: %d numTrips: %d\n", length, numTriplets);
	// allocate the triplets adding one to account for the null byte
	m_pSeq = new char[numBytes];
	memset(m_pSeq, 0, numBytes);
	
	const char* strData = seq.data();
	
	for(int i = 0; i < m_length; i++)
	{
		int byteNumber = seqIndexToByteNumber(i);
		int baseIndex = seqIndexToBaseIndex(i);
		
		//printf("setting %d = (%d %d) to %c\n", i, byteNumber, baseIndex, strData[i]);
		setBase(m_pSeq, byteNumber, baseIndex, strData[i]);
	}
}

// Copy constructor
PackedSeq::PackedSeq(const PackedSeq& pseq)
{
	// allocate memory and copy over the seq
	m_length = pseq.m_length;
	int numBytes = getNumCodingBytes(m_length);
	m_pSeq = new char[numBytes];
	memset(m_pSeq, 0, numBytes);	
	
	// copy the sequence over
	memcpy(m_pSeq, pseq.m_pSeq, numBytes);
}

// Assignment operator
PackedSeq& PackedSeq::operator=(const PackedSeq& other)
{
	// Detect self assignment
	if (this == &other)
	{
		return *this;
	}
	
	// Delete previous seq, this will either be NULL or valid allocated memory
	delete [] m_pSeq;
	m_pSeq = 0;
	
	// Allocate and fill the new seq
	m_length = other.m_length;	
	int numBytes = getNumCodingBytes(m_length);
	
	m_pSeq = new char[numBytes];
	memset(m_pSeq, 0, numBytes);	
	
	// copy the sequence over
	memcpy(m_pSeq, other.m_pSeq, numBytes);
	
	return *this;
}

bool PackedSeq::operator==(const PackedSeq& other) const
{
	int numBytes = getNumCodingBytes(m_length);
	return (memcmp(m_pSeq, other.m_pSeq, numBytes) == 0);
}

bool PackedSeq::operator!=(const PackedSeq& other) const
{
	return !(*this == other);
}

bool PackedSeq::operator<(const PackedSeq& other) const
{
	int numBytes = getNumCodingBytes(m_length);
	return (memcmp(m_pSeq, other.m_pSeq, numBytes) < 0);
}

char PackedSeq::getBase(int seqIndex) const
{
	assert(seqIndex < m_length);
	int byteNumber = seqIndexToByteNumber(seqIndex);
	int baseIndex = seqIndexToBaseIndex(seqIndex);	
	return getBase(m_pSeq, byteNumber, baseIndex);
}


int PackedSeq::getSequenceLength() const
{
	return m_length;
}

const char* const PackedSeq::getDataPtr() const
{
	return m_pSeq;	
}

void PackedSeq::print() const
{
	const char* iter = m_pSeq;
	while(*iter)
	{
		printf("%c", *iter);
		iter++;
	}
}

// Get the number of coding triplets
int PackedSeq::getNumCodingBytes(int seqLength)
{
	if(seqLength % 4 == 0)
	{
		return seqLength / 4;
	}
	else
	{
		return (seqLength / 4) + 1;
	} 
}

Sequence PackedSeq::decode() const
{
	// allocate space for the new string
	Sequence outstr;
	
	for(int i = 0; i < m_length; i++)
	{

		char base = getBase(i);
		//printf("decoding (%d %d) to %c\n", tripletNumber, baseIndex, base);
		outstr.push_back(base);
	}
	
	//printf("decode: %s\n", outstr.c_str());
	return outstr; 
}

// Change this sequence into its reverse complement

//TODO: further optimize this. 1) don't need to go code->base->code 2) can reverse within bytes then swap whole bytes?
void PackedSeq::reverseComplement()
{
	int numBytes = getNumCodingBytes(m_length);

	// reverse the string by swapping
	for(int i = 0; i < m_length / 2; i++)
	{
		int revPos = m_length - i - 1;
		
		int read1ByteNumber = seqIndexToByteNumber(i);
		int read1BaseIndex = seqIndexToBaseIndex(i);
		
		int read2ByteNumber = seqIndexToByteNumber(revPos);
		int read2BaseIndex = seqIndexToBaseIndex(revPos);		
		
		char base1 = getBase(m_pSeq, read1ByteNumber, read1BaseIndex);
		char base2 = getBase(m_pSeq, read2ByteNumber, read2BaseIndex);
		setBase(m_pSeq, read1ByteNumber, read1BaseIndex, base2);
		setBase(m_pSeq, read2ByteNumber, read2BaseIndex, base1);
	}
	
	for(int i = 0; i < numBytes; i++)
	{
		// complement each byte
		m_pSeq[i] = ~m_pSeq[i];
	}
}

char PackedSeq::shiftAppend(char base)
{
	// shift the sequence left and append a new base to the end
	int numBytes = getNumCodingBytes(m_length);

	char shiftIn = base;
	
	// starting from the last byte, shift the new base in and get the captured base
	for(int i = numBytes - 1; i >= 0; i--)
	{
		// calculate the index
		// if this is the last byte, use 
		int index = (i == (numBytes - 1)) ? seqIndexToBaseIndex(m_length - 1) : 3;
		shiftIn = leftShiftByte(m_pSeq, i, index, shiftIn);
	}
	
	// return the base shifted out of the first byte
	return shiftIn;
}

char PackedSeq::shiftPrepend(char base)
{
	// shift the sequence right and append a new base to the end
	int numBytes = getNumCodingBytes(m_length);

	int lastBaseByte = seqIndexToByteNumber(m_length - 1);
	int lastBaseIndex = seqIndexToBaseIndex(m_length - 1);
	
	// save the last base (which gets shifted out)
	char lastBase = getBase(m_pSeq, lastBaseByte, lastBaseIndex);
	
	char shiftIn = base;
	
	// starting from the last byte, shift the new base in and get the captured base
	for(int i = 0; i <= numBytes - 1; i++)
	{
		// index is always zero
		int index = 0;
		shiftIn = rightShiftByte(m_pSeq, i, index, shiftIn);
	}
		
	return lastBase;	
}

char PackedSeq::leftShiftByte(char* pSeq, int byteNum, int index, char base)
{
	// save the first base
	char outBase = (pSeq[byteNum] >> 6) & 0x3;
	
	// shift left one position
	pSeq[byteNum] <<= 2;
	
	// Set the new base
	setBase(pSeq, byteNum, index, base);

	return codeToBase(outBase);
}

char PackedSeq::rightShiftByte(char* pSeq, int byteNum, int index, char base)
{
	// save the last base
	char outBase = pSeq[byteNum] & 0x3;
	
	// shift right one position
	pSeq[byteNum] >>= 2;
	
	// add the new base
	setBase(pSeq, byteNum, index, base);
	
	return codeToBase(outBase);
}

// Set/Get a particular base
void PackedSeq::setBase(char* pSeq, int byteNum, int index, char base)
{
	//printf("setting (%d %d) to %c\n", byteNum, index, base);
	char val = baseToCode(base);
	
	// shift the value into position
	int shiftValue = 2*(3 - index);
	val <<= shiftValue;
	
	// clear the value
	char mask = 0x3;
	mask <<= shiftValue;
	mask = ~mask;
	pSeq[byteNum] &= mask;
	
	
	// set the appropriate value with an OR
	pSeq[byteNum] |= val;
}

char PackedSeq::getBase(const char* pSeq, int byteNum, int index) const
{
	int shiftLen = 2 * (3 - index);
	return codeToBase((pSeq[byteNum] >> shiftLen) & 0x3);
}

int PackedSeq::seqIndexToByteNumber(int seqIndex) const
{
	return seqIndex / 4;
}

int PackedSeq::seqIndexToBaseIndex(int seqIndex) const
{
	return seqIndex % 4; 
}

// return the two bit code for each base
// the input base MUST be in upper case
char PackedSeq::baseToCode(char base) const
{
	if(base == 'A')
	{
		return 0;
	}
	else if(base == 'C')
	{
		return 1;
	}
	else if(base == 'G')
	{
		return 2;
	}
	else if(base == 'T')
	{
		return 3;
	}
	else
	{
		// unknown base
		assert(false);
		return 0;
	}
}

char PackedSeq::codeToBase(char code) const
{
	if(code == 0)
	{
		return 'A';
	}
	else if(code == 1)
	{
		return 'C';
	}
	else if(code == 2)
	{
		return 'G';
	}
	else if(code == 3)
	{
		return 'T';
	}
	else
	{
		// unknown code
		assert(false);
		return 'A';
	}
}

PackedSeq reverseComplement(const PackedSeq& seq)
{
	PackedSeq rc(seq);
	rc.reverseComplement();
	return rc;	
}
