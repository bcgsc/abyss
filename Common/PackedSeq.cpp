#include <math.h>
#include "PackedSeq.h"

//
// Default constructor
//
PackedSeq::PackedSeq() : m_length(0), m_flags(0)
{
	memset(m_seq, 0, NUM_BYTES);
}

//
// Destructor
//
PackedSeq::~PackedSeq()
{

}

//
// Construct a sequence from a series of bytes
//
PackedSeq::PackedSeq(const char* const pData, int length) : m_flags(0)
{
	assert(length < MAX_KMER);
	memset(m_seq, 0, NUM_BYTES);
	int numBytes = getNumCodingBytes(length);

	// copy over the bytes
	memcpy(m_seq, pData, numBytes);
	
	// set the sequence length
	m_length = length;
}

//
// Construct a sequence from a String-based sequence
//
PackedSeq::PackedSeq(const Sequence& seq) : m_flags(0)
{
	memset(m_seq, 0, NUM_BYTES);
	int length = seq.length();
	
	assert(length < MAX_KMER);
	
	//printf("packing %s\n", seq.c_str());
	// calculate the number of triplets required

	m_length = length;

	const char* strData = seq.data();
	
	for(int i = 0; i < m_length; i++)
	{
		setBase(m_seq, i, strData[i]);
	}
}

//
// Copy constructor
//
PackedSeq::PackedSeq(const PackedSeq& pseq)
{
	memset(m_seq, 0, NUM_BYTES);
	// allocate memory and copy over the seq
	m_length = pseq.m_length;
	int numBytes = getNumCodingBytes(m_length);
	
	// copy the sequence over
	memcpy(m_seq, pseq.m_seq, numBytes);
	
	m_flags = pseq.m_flags;
	m_extensions[0] = pseq.m_extensions[0];
	m_extensions[1] = pseq.m_extensions[1];
}

//
// Assignment operator
//
PackedSeq& PackedSeq::operator=(const PackedSeq& other)
{
	// Detect self assignment
	if (this == &other)
	{
		return *this;
	}
	memset(m_seq, 0, NUM_BYTES);
	// Allocate and fill the new seq
	m_length = other.m_length;	
	int numBytes = getNumCodingBytes(m_length);

	// copy the sequence over
	memcpy(m_seq, other.m_seq, numBytes);
	
	m_flags = other.m_flags;
	
	return *this;
}

//
// Equality operator
//
bool PackedSeq::operator==(const PackedSeq& other) const
{
	return compare(other) == 0;
}

int PackedSeq::compare(const PackedSeq& other) const
{
	assert(m_length == other.m_length);
	// This is valid since the remaining portion of the buffer is memset to zero
	// This HAS to be the case or else the equality will not work
	//return memcmp(m_seq, other.m_seq, NUM_BYTES);
	
	int numBytes = getNumCodingBytes(m_length);
	if(m_length % 4 != 0)
	{
		numBytes--;
	}
	
	int firstCompare = memcmp(m_seq, other.m_seq, numBytes);
	
	if(firstCompare != 0)
	{
		return firstCompare;
	}
	else
	{
		for(int i = 4*numBytes; i < m_length; i++)
		{
			
			char base1 = this->getBase(i);
			char base2 = other.getBase(i);
			if(base1 != base2)
			{
				return base1 - base2;
			}
		}
	}
	return 0;
	
}

//
// Inequality operator
//
bool PackedSeq::operator!=(const PackedSeq& other) const
{
	return !(*this == other);
}

//
// Less than operator for the strings, performs an O(n) comparison (TODO: can be optimized??)
// this is slow as hell, fix it
//
bool PackedSeq::operator<(const PackedSeq& other) const
{	
	return compare(other) < 0;
}

//
// return a subsequence of this sequence
//
PackedSeq PackedSeq::subseq(int start, int len) const
{
	// allocate space to hold the temporary string
	int numBytes = getNumCodingBytes(len);
	char* tempBuffer = new char[numBytes];
	memset(tempBuffer, 0, numBytes);
	
	for(int i = start; i < start + len; i++)
	{
		int index = i - start;
		setBase(tempBuffer, index, getBase(i));
	}
	
	PackedSeq sub(tempBuffer, len);
	delete [] tempBuffer;
	return sub;
}

//
// Get the length of the sequence
//
int PackedSeq::getSequenceLength() const
{
	return m_length;
}

//
// Return a pointer to the raw data
//
const char* const PackedSeq::getDataPtr() const
{
	return m_seq;	
}

//
// Print the string to stdout
//
void PackedSeq::print() const
{
	const char* iter = getDataPtr();
	while(*iter)
	{
		printf("%c", *iter);
		iter++;
	}
}

//
// Get the number of coding bytes
// 
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

//
// This function computes a hash-like value of the packed sequence over the first 8 bases and the reverse complement of the last 8 bases
// The reverse complement of the last 8 bases is used so that a sequence and its reverse-comp will hash to the same value which is desirable in this case
// Todo: make this faster?
//
unsigned int PackedSeq::getCode() const
{
	const int NUM_BYTES = 4;
	char firstByte[NUM_BYTES];
	for(int i  = 0; i < NUM_BYTES; i++)
	{
		firstByte[i] = m_seq[i];
	}
	
	// Fill out the last byte as the reverse complement of the last 8 bases
	char lastByte[NUM_BYTES];
	int prime = 101;
	
	for(int i = 0; i < 4*NUM_BYTES; i++)
	{
		int index = m_length - 1 - i;
		char base = getBase(index);
		setBase(lastByte, i/4, i % 4, complement(base));
	}

	char f[NUM_BYTES];
	unsigned int sum = 0;
	
	for(int i = 0; i < NUM_BYTES; i++)
	{
		int total = 1;
		for(int j = 0; j < i; j++)
		{
			total *= prime;
		}
		
		sum += (firstByte[i] ^ lastByte[i]) * total;
	}
	
	//short f1 = (firstByte[1] ^ lastByte[1]);
	//printf("%d code\n", code);
	return sum;
}

//
// Decode this sequence into an ascii string
//
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

//
// Decode a single byte into it's ascii representation
//
Sequence PackedSeq::decodeByte(const char byte) const
{
	Sequence outstr;
	outstr.reserve(4);
	for(int i = 0; i < 4; i++)
	{
		char base = getBase(&byte, 0, i);
		outstr.push_back(base);
	}
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
		
		char base1 = getBase(i);
		char base2 = getBase(revPos);
		setBase(m_seq, i, base2);
		setBase(m_seq, revPos, base1);
	}
	
	for(int i = 0; i < numBytes; i++)
	{
		// complement each byte
		m_seq[i] = ~m_seq[i];
	}
}

//
//
//
char PackedSeq::rotate(extDirection dir, char base)
{
	if(dir == SENSE)
	{
		return shiftAppend(base);
	}
	else
	{
		return shiftPrepend(base);	
	}
}

//
//
//
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
		shiftIn = leftShiftByte(m_seq, i, index, shiftIn);
	}
	
	// return the base shifted out of the first byte
	return shiftIn;
}

//
//
//
char PackedSeq::shiftPrepend(char base)
{
	// shift the sequence right and append a new base to the end
	int numBytes = getNumCodingBytes(m_length);

	int lastBaseByte = seqIndexToByteNumber(m_length - 1);
	int lastBaseIndex = seqIndexToBaseIndex(m_length - 1);
	
	// save the last base (which gets shifted out)
	char lastBase = getBase(m_seq, lastBaseByte, lastBaseIndex);
	
	char shiftIn = base;
	
	// starting from the last byte, shift the new base in and get the captured base
	for(int i = 0; i <= numBytes - 1; i++)
	{
		// index is always zero
		int index = 0;
		shiftIn = rightShiftByte(m_seq, i, index, shiftIn);
	}
		
	return lastBase;	
}

//
//
//
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

//
//
//
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

//
//
//
void PackedSeq::setExtension(extDirection dir, char b)
{
	m_extensions[dir].SetBase(b);
}


//
//
//
void PackedSeq::clearExtension(extDirection dir, char b)
{
	m_extensions[dir].ClearBase(b);
}

//
//
//
bool PackedSeq::checkExtension(extDirection dir, char b) const
{
	return m_extensions[dir].CheckBase(b);	
}

//
//
//
bool PackedSeq::hasExtension(extDirection dir) const
{
	return m_extensions[dir].HasExtension();	
}

//
//
//
bool PackedSeq::isAmbiguous(extDirection dir) const
{
	return m_extensions[dir].IsAmbiguous();	
}

//
//
//
void PackedSeq::printExtension() const
{
	printf("seq: %s\n", decode().c_str());
	printf("sxt: ");
	m_extensions[SENSE].print();
	
	printf("as : ");
	m_extensions[ANTISENSE].print();	
		
}

//
//
//
void PackedSeq::setFlag(SeqFlag flag)
{
	m_flags |= flag;		
}

//
//
//
bool PackedSeq::isFlagSet(SeqFlag flag) const
{
	return m_flags & flag;
}

//
// set a base by the actual index [0, length)
// beware, this does not check for out of bounds access
//
void PackedSeq::setBase(char* pSeq, int seqIndex, char base)
{
	int byteNumber = seqIndexToByteNumber(seqIndex);
	int baseIndex = seqIndexToBaseIndex(seqIndex);	
	return setBase(pSeq, byteNumber, baseIndex, base);
}

//
//Set a base by byte number/ sub index
// beware, this does not check for out of bounds access
//
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

//
// get a base by the actual index [0, length)
//
char PackedSeq::getBase(int seqIndex) const
{
	assert(seqIndex < m_length);
	int byteNumber = seqIndexToByteNumber(seqIndex);
	int baseIndex = seqIndexToBaseIndex(seqIndex);	
	return getBase(m_seq, byteNumber, baseIndex);
}

//
// get a base by the byte number and sub index
//
char PackedSeq::getBase(const char* pSeq, int byteNum, int index) const
{
	int shiftLen = 2 * (3 - index);
	return codeToBase((pSeq[byteNum] >> shiftLen) & 0x3);
}

//
//
//
int PackedSeq::seqIndexToByteNumber(int seqIndex)
{
	return seqIndex / 4;
}

//
//
//
int PackedSeq::seqIndexToBaseIndex(int seqIndex)
{
	return seqIndex % 4; 
}

//
// return the two bit code for each base
// the input base MUST be in upper case
//

char PackedSeq::baseToCode(char base)
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

//
//
//
char PackedSeq::codeToBase(char code)
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

//
//
//
PackedSeq reverseComplement(const PackedSeq& seq)
{
	PackedSeq rc(seq);
	rc.reverseComplement();
	return rc;	
}
