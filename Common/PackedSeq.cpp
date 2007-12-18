#include <math.h>
#include "PackedSeq.h"

// Construct a sequence with the memory already allocated
// The class is responsible for freeing the data
PackedSeq::PackedSeq(char* const pData, int length)
{
	m_length = length;
	m_pSeq = pData;
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
		setBase(m_pSeq, byteNumber, baseIndex, strData[i]);
	}
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
	
	//printf("num total bases: %d\n", numBasesTotal);
	// allocate space for the new string
	Sequence outstr;
	
	for(int i = 0; i < m_length; i++)
	{
		int byteNumber = seqIndexToByteNumber(i);
		int baseIndex = seqIndexToBaseIndex(i);
		char base = getBase(m_pSeq, byteNumber, baseIndex);
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

PackedSeq::~PackedSeq()
{
	delete [] m_pSeq;
	m_pSeq = NULL;	
}
