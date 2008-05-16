#include <math.h>
#include "PackedSeq.h"

//
// Default constructor
//
PackedSeq::PackedSeq() : m_length(0), m_flags(0), m_multiplicity(1)
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
PackedSeq::PackedSeq(const char* const pData, int length) : m_flags(0), m_multiplicity(1)
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
PackedSeq::PackedSeq(const Sequence& seq) : m_flags(0), m_multiplicity(1)
{
	memset(m_seq, 0, NUM_BYTES);
	int length = seq.length();

	assert(length <= MAX_KMER);
	
	//printf("packing %s\n", seq.c_str());
	// calculate the number of triplets required

	m_length = length;

	const char* strData = seq.data();
	
	for(int i = 0; i < m_length; i++)
		setBaseChar(m_seq, i, strData[i]);
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
	m_extRecord.dir[SENSE] = pseq.m_extRecord.dir[SENSE];
	m_extRecord.dir[ANTISENSE] = pseq.m_extRecord.dir[ANTISENSE];
	m_multiplicity = pseq.m_multiplicity;
}

//
// Serialize this packed seq to the buffer
// TODO: make this machine independent (using mpi datatypes?)
// Return the number of bytes read
//
size_t PackedSeq::serialize(char* buffer) const
{
	size_t offset = 0;
	
	memcpy(buffer + offset, &m_seq, sizeof(m_seq));
	offset += sizeof(m_seq);
	
	memcpy(buffer + offset, &m_length, sizeof(m_length));
	offset += sizeof(m_length);
	
	memcpy(buffer + offset, &m_flags, sizeof(m_flags));
	offset += sizeof(m_flags);
	
	memcpy(buffer + offset, &m_multiplicity, sizeof(m_multiplicity));
	offset += sizeof(m_multiplicity);	
	
	memcpy(buffer + offset, &m_extRecord, sizeof(m_extRecord));
	offset += sizeof(m_extRecord);
	
	assert(offset == sizeof(PackedSeq));	

	return offset;		
}

//
// Unserialize this packed seq from the buffer
// TODO: make this machine independent (using mpi datatypes?)
//
size_t PackedSeq::unserialize(const char* buffer)
{
	size_t offset = 0;
	
	memcpy(m_seq, buffer + offset, sizeof(m_seq));
	offset += sizeof(m_seq);
	
	memcpy(&m_length, buffer + offset, sizeof(m_length));
	offset += sizeof(m_length);
	
	memcpy(&m_flags, buffer + offset, sizeof(m_flags));
	offset += sizeof(m_flags);
	
	memcpy(&m_multiplicity, buffer + offset, sizeof(m_multiplicity));
	offset += sizeof(m_multiplicity);	
	
	memcpy(&m_extRecord, buffer + offset, sizeof(m_extRecord));
	offset += sizeof(m_extRecord);
	
	assert(offset == sizeof(PackedSeq));
	
	return offset;			
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
	
	
	m_extRecord.dir[SENSE] = other.m_extRecord.dir[SENSE];
	m_extRecord.dir[ANTISENSE] = other.m_extRecord.dir[ANTISENSE];
	m_multiplicity = other.m_multiplicity;

	return *this;
}

//
// Compare two sequences.
// 
int PackedSeq::compare(const PackedSeq& other) const
{
	assert(m_length == other.m_length);
	unsigned nbytes = getNumCodingBytes(m_length);
	return memcmp(m_seq, other.m_seq, nbytes);
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
unsigned PackedSeq::getNumCodingBytes(unsigned seqLength)
{
	return (seqLength + 3) / 4;
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
		setBaseCode(lastByte, i/4, i % 4,
				complementBaseCode(getBaseCode(index)));
	}

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

uint32_t hashlittle( const void *key, size_t length, uint32_t initval);	
size_t PackedSeq::getHashCode() const
{
	// Hash on the numbytes - 1. This is to avoid getting different hash values for the same sequence for n % 4 != 0 sequences
	int code = hashlittle(m_seq, getNumCodingBytes(m_length) - 1, 131);

	return code;
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

		char base = getBaseChar(i);
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
		char base = getBaseChar(&byte, 0, i);
		outstr.push_back(base);
	}
	return outstr;
}

/** Swap the positions of four bases. */
static const uint8_t swapBases[256] = {
	0x00, 0x40, 0x80, 0xc0, 0x10, 0x50, 0x90, 0xd0,
	0x20, 0x60, 0xa0, 0xe0, 0x30, 0x70, 0xb0, 0xf0,
	0x04, 0x44, 0x84, 0xc4, 0x14, 0x54, 0x94, 0xd4,
	0x24, 0x64, 0xa4, 0xe4, 0x34, 0x74, 0xb4, 0xf4,
	0x08, 0x48, 0x88, 0xc8, 0x18, 0x58, 0x98, 0xd8,
	0x28, 0x68, 0xa8, 0xe8, 0x38, 0x78, 0xb8, 0xf8,
	0x0c, 0x4c, 0x8c, 0xcc, 0x1c, 0x5c, 0x9c, 0xdc,
	0x2c, 0x6c, 0xac, 0xec, 0x3c, 0x7c, 0xbc, 0xfc,
	0x01, 0x41, 0x81, 0xc1, 0x11, 0x51, 0x91, 0xd1,
	0x21, 0x61, 0xa1, 0xe1, 0x31, 0x71, 0xb1, 0xf1,
	0x05, 0x45, 0x85, 0xc5, 0x15, 0x55, 0x95, 0xd5,
	0x25, 0x65, 0xa5, 0xe5, 0x35, 0x75, 0xb5, 0xf5,
	0x09, 0x49, 0x89, 0xc9, 0x19, 0x59, 0x99, 0xd9,
	0x29, 0x69, 0xa9, 0xe9, 0x39, 0x79, 0xb9, 0xf9,
	0x0d, 0x4d, 0x8d, 0xcd, 0x1d, 0x5d, 0x9d, 0xdd,
	0x2d, 0x6d, 0xad, 0xed, 0x3d, 0x7d, 0xbd, 0xfd,
	0x02, 0x42, 0x82, 0xc2, 0x12, 0x52, 0x92, 0xd2,
	0x22, 0x62, 0xa2, 0xe2, 0x32, 0x72, 0xb2, 0xf2,
	0x06, 0x46, 0x86, 0xc6, 0x16, 0x56, 0x96, 0xd6,
	0x26, 0x66, 0xa6, 0xe6, 0x36, 0x76, 0xb6, 0xf6,
	0x0a, 0x4a, 0x8a, 0xca, 0x1a, 0x5a, 0x9a, 0xda,
	0x2a, 0x6a, 0xaa, 0xea, 0x3a, 0x7a, 0xba, 0xfa,
	0x0e, 0x4e, 0x8e, 0xce, 0x1e, 0x5e, 0x9e, 0xde,
	0x2e, 0x6e, 0xae, 0xee, 0x3e, 0x7e, 0xbe, 0xfe,
	0x03, 0x43, 0x83, 0xc3, 0x13, 0x53, 0x93, 0xd3,
	0x23, 0x63, 0xa3, 0xe3, 0x33, 0x73, 0xb3, 0xf3,
	0x07, 0x47, 0x87, 0xc7, 0x17, 0x57, 0x97, 0xd7,
	0x27, 0x67, 0xa7, 0xe7, 0x37, 0x77, 0xb7, 0xf7,
	0x0b, 0x4b, 0x8b, 0xcb, 0x1b, 0x5b, 0x9b, 0xdb,
	0x2b, 0x6b, 0xab, 0xeb, 0x3b, 0x7b, 0xbb, 0xfb,
	0x0f, 0x4f, 0x8f, 0xcf, 0x1f, 0x5f, 0x9f, 0xdf,
	0x2f, 0x6f, 0xaf, 0xef, 0x3f, 0x7f, 0xbf, 0xff
};

struct Seq {
	uint64_t x[(PackedSeq::NUM_BYTES + 7)/8];
};

/** Load with appropriate endianness for shifting. */
static Seq load(const uint8_t *src)
{
	assert(PackedSeq::NUM_BYTES == 10);
	Seq seq = {
		(uint64_t)src[0] << 56
			| (uint64_t)src[1] << 48
			| (uint64_t)src[2] << 40
			| (uint64_t)src[3] << 32
			| (uint64_t)src[4] << 24
			| (uint64_t)src[5] << 16
			| (uint64_t)src[6] << 8
			| (uint64_t)src[7] << 0,
		(uint64_t)src[8] << 56
			| (uint64_t)src[9] << 48
	};
	return seq;
}

static void store(uint8_t *dest, Seq seq)
{
	dest[0] = seq.x[0] >> 56;
	dest[1] = seq.x[0] >> 48;
	dest[2] = seq.x[0] >> 40;
	dest[3] = seq.x[0] >> 32;
	dest[4] = seq.x[0] >> 24;
	dest[5] = seq.x[0] >> 16;
	dest[6] = seq.x[0] >> 8;
	dest[7] = seq.x[0] >> 0;
	dest[8] = seq.x[1] >> 56;
	dest[9] = seq.x[1] >> 48;
}

/**
 * Reverse the bytes by storing them in the reverse order of
 * loading, and reverse the words in the same fasion.
 */
static void storeReverse(uint8_t *dest, Seq seq)
{
	uint64_t *p0 = (uint64_t *)&dest[0];
	uint16_t *p1 = (uint16_t *)&dest[8];
	*p0 = seq.x[1];
	*p1 = seq.x[0];
}

/** Shift right by the specified number of bits. */
static void shiftRight(Seq *pseq, uint8_t n)
{
	if (n == 0)
		return;
	uint64_t x0 = pseq->x[0], x1 = pseq->x[1];
	if (n < 64) {
		pseq->x[0] = x0 >> n;
		pseq->x[1] = x1 >> n | x0 << (64 - n);
	} else {
		pseq->x[0] = 0;
		pseq->x[1] = x0 >> (n - 64);
	}
}

/** Shift left by the specified number of bits. */
static void shiftLeft(Seq *pseq, uint8_t n)
{
	if (n == 0)
		return;
	uint64_t x0 = pseq->x[0], x1 = pseq->x[1];
	if (n < 64) {
		pseq->x[0] = x0 << n | x1 >> (64 - n);
		pseq->x[1] = x1 << n;
	} else {
		pseq->x[0] = x1 << (n - 64);
		pseq->x[1] = 0;
	}
}

//
// Return a subsequence of this sequence.
//
PackedSeq PackedSeq::subseq(int start, int len) const
{
	assert(start >= 0 && len >= 0);
	assert(start + len <= m_length);
	Seq seq = load((uint8_t*)m_seq);
#if 0
	// A simple left shift would suffice if it weren't necessary to
	// zero the bases beyond the end of the sequence. compare requires
	// that they are zeroed.
	shiftLeft(&seq, 2*start);
#else
	// Zero the bases beyond the end of the sequence. A mask would
	// probably be faster.
	shiftRight(&seq, 128 - 2*(len + start));
	shiftLeft(&seq, 128 - 2*len);
#endif
	PackedSeq sub;
	sub.m_length = len;
	store((uint8_t*)sub.m_seq, seq);
	return sub;
}

//
// Change this sequence into its reverse complement
//
void PackedSeq::reverseComplement()
{
	Seq seq = load((uint8_t*)m_seq);

	// Complement the bits.
	seq.x[0] = ~seq.x[0];
	seq.x[1] = ~seq.x[1];

	// Shift the bits flush to the right of the double word.
	shiftRight(&seq, 128 - 2*m_length);

	storeReverse((uint8_t*)m_seq, seq);

	// Reverse the pairs of bits withing a byte.
	unsigned numBytes = getNumCodingBytes(m_length);
	for (unsigned i = 0; i < numBytes; i++)
		m_seq[i] = swapBases[(uint8_t)m_seq[i]];
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

void PackedSeq::setLastBase(extDirection dir, char base)
{
	setBaseChar(m_seq, dir == SENSE ? m_length - 1 : 0, base);
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
	char lastBase = getBaseChar(m_seq, lastBaseByte, lastBaseIndex);
	
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
	setBaseChar(pSeq, byteNum, index, base);

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
	setBaseChar(pSeq, byteNum, index, base);
	
	return codeToBase(outBase);
}

//
//
//
void PackedSeq::setExtension(extDirection dir, SeqExt extension)
{
	m_extRecord.dir[dir] = extension;
}

//
//
//
void PackedSeq::setBaseExtension(extDirection dir, char b)
{
	m_extRecord.dir[dir].SetBase(b);
}

//
//
//
void PackedSeq::clearAllExtensions(extDirection dir)
{
	m_extRecord.dir[dir].ClearAll();
}


//
//
//
void PackedSeq::clearExtension(extDirection dir, char b)
{
	m_extRecord.dir[dir].ClearBase(b);
}

//
//
//
bool PackedSeq::checkExtension(extDirection dir, char b) const
{
	return m_extRecord.dir[dir].CheckBase(b);	
}

//
//
//
bool PackedSeq::hasExtension(extDirection dir) const
{
	return m_extRecord.dir[dir].HasExtension();	
}

//
//
//
bool PackedSeq::isAmbiguous(extDirection dir) const
{
	return m_extRecord.dir[dir].IsAmbiguous();	
}

//
// Return the sequences extension in the specified direction
//
SeqExt PackedSeq::getExtension(extDirection dir) const
{
	return m_extRecord.dir[dir];
}

//
//
//
void PackedSeq::printExtension() const
{
	printf("seq: %s\n", decode().c_str());
	printf("sxt: ");
	m_extRecord.dir[SENSE].print();
	
	printf("as : ");
	m_extRecord.dir[ANTISENSE].print();	
		
}

//
// set a base by the index [0, length)
// beware, this does not check for out of bounds access
//
void PackedSeq::setBaseCode(char* pSeq, int seqIndex, uint8_t base)
{
	int byteNumber = seqIndexToByteNumber(seqIndex);
	int baseIndex = seqIndexToBaseIndex(seqIndex);	
	return setBaseCode(pSeq, byteNumber, baseIndex, base);
}

//
// set a base by the index [0, length)
// beware, this does not check for out of bounds access
//
void PackedSeq::setBaseChar(char* pSeq, int seqIndex, char base)
{
	return setBaseCode(pSeq, seqIndex, baseToCode(base));
}

//
//Set a base by byte number/ sub index
// beware, this does not check for out of bounds access
//
void PackedSeq::setBaseCode(char* pSeq, int byteNum, int index,
		uint8_t base)
{
	// shift the value into position
	int shiftValue = 2*(3 - index);
	base <<= shiftValue;
	
	// clear the value
	char mask = 0x3;
	mask <<= shiftValue;
	mask = ~mask;
	pSeq[byteNum] &= mask;
	
	
	// set the appropriate value with an OR
	pSeq[byteNum] |= base;
}

//
//Set a base by byte number/ sub index
// beware, this does not check for out of bounds access
//
void PackedSeq::setBaseChar(char* pSeq, int byteNum, int index, char base)
{
	setBaseCode(pSeq, byteNum, index, baseToCode(base));
}

//
// get a base code by the index [0, length)
//
uint8_t PackedSeq::getBaseCode(int seqIndex) const
{
	assert(seqIndex < m_length);
	int byteNumber = seqIndexToByteNumber(seqIndex);
	int baseIndex = seqIndexToBaseIndex(seqIndex);	
	return getBaseCode(m_seq, byteNumber, baseIndex);
}

//
// get a base by the index [0, length)
//
char PackedSeq::getBaseChar(int seqIndex) const
{
	return codeToBase(getBaseCode(seqIndex));
}

//
// get a base code by the byte number and sub index
//
uint8_t PackedSeq::getBaseCode(const char* pSeq, int byteNum, int index) const
{
	int shiftLen = 2 * (3 - index);
	return (pSeq[byteNum] >> shiftLen) & 0x3;
}

//
// get a base by the byte number and sub index
//
char PackedSeq::getBaseChar(const char* pSeq, int byteNum, int index) const
{
	return codeToBase(getBaseCode(pSeq, byteNum, index));
}

//
//
//
unsigned PackedSeq::seqIndexToByteNumber(unsigned seqIndex)
{
	return seqIndex / 4;
}

//
//
//
unsigned PackedSeq::seqIndexToBaseIndex(unsigned seqIndex)
{
	return seqIndex % 4; 
}

//
// return the two bit code for each base
// the input base MUST be in upper case
//
uint8_t PackedSeq::baseToCode(char base)
{
	unsigned i = base - 'A';
	assert(i < 26);
	static const int8_t table[26] = {
		 0, -1,  1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1
	};
	assert(table[i] >= 0);
	return table[i];
}

//
//
//
char PackedSeq::codeToBase(uint8_t code)
{
	assert(code < 4);
	return "ACGT"[code];
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
/*
-------------------------------------------------------------------------------
lookup3.c, by Bob Jenkins, May 2006, Public Domain.

These are functions for producing 32-bit hashes for hash table lookup.
hashword(), hashlittle(), hashlittle2(), hashbig(), mix(), and final() 
are externally useful functions.  Routines to test the hash are included 
if SELF_TEST is defined.  You can use this free for any purpose.  It's in
the public domain.  It has no warranty.

You probably want to use hashlittle().  hashlittle() and hashbig()
hash byte arrays.  hashlittle() is is faster than hashbig() on
little-endian machines.  Intel and AMD are little-endian machines.
On second thought, you probably want hashlittle2(), which is identical to
hashlittle() except it returns two 32-bit hashes for the price of one.  
You could implement hashbig2() if you wanted but I haven't bothered here.

If you want to find a hash of, say, exactly 7 integers, do
  a = i1;  b = i2;  c = i3;
  mix(a,b,c);
  a += i4; b += i5; c += i6;
  mix(a,b,c);
  a += i7;
  final(a,b,c);
then use c as the hash value.  If you have a variable length array of
4-byte integers to hash, use hashword().  If you have a byte array (like
a character string), use hashlittle().  If you have several byte arrays, or
a mix of things, see the comments above hashlittle().  

Why is this so big?  I read 12 bytes at a time into 3 4-byte integers, 
then mix those integers.  This is fast (you can do a lot more thorough
mixing with 12*3 instructions on 3 integers than you can with 3 instructions
on 1 byte), but shoehorning those bytes into integers efficiently is messy.
-------------------------------------------------------------------------------
*/
#define SELF_TEST 1

#include <stdio.h>      /* defines printf for tests */
#include <time.h>       /* defines time_t for timings in the test */
#include <stdint.h>     /* defines uint32_t etc */
#include <sys/param.h>  /* attempt to define endianness */
#ifdef linux
# include <endian.h>    /* attempt to define endianness */
#endif

/*
 * My best guess at if you are big-endian or little-endian.  This may
 * need adjustment.
 */
#if (defined(__BYTE_ORDER) && defined(__LITTLE_ENDIAN) && \
     __BYTE_ORDER == __LITTLE_ENDIAN) || \
    (defined(i386) || defined(__i386__) || defined(__i486__) || \
     defined(__i586__) || defined(__i686__) || defined(vax) || defined(MIPSEL))
# define HASH_LITTLE_ENDIAN 1
# define HASH_BIG_ENDIAN 0
#elif (defined(__BYTE_ORDER) && defined(__BIG_ENDIAN) && \
       __BYTE_ORDER == __BIG_ENDIAN) || \
      (defined(sparc) || defined(POWERPC) || defined(mc68000) || defined(sel))
# define HASH_LITTLE_ENDIAN 0
# define HASH_BIG_ENDIAN 1
#else
# define HASH_LITTLE_ENDIAN 0
# define HASH_BIG_ENDIAN 0
#endif

#define hashsize(n) ((uint32_t)1<<(n))
#define hashmask(n) (hashsize(n)-1)
#define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))

/*
-------------------------------------------------------------------------------
mix -- mix 3 32-bit values reversibly.

This is reversible, so any information in (a,b,c) before mix() is
still in (a,b,c) after mix().

If four pairs of (a,b,c) inputs are run through mix(), or through
mix() in reverse, there are at least 32 bits of the output that
are sometimes the same for one pair and different for another pair.
This was tested for:
* pairs that differed by one bit, by two bits, in any combination
  of top bits of (a,b,c), or in any combination of bottom bits of
  (a,b,c).
* "differ" is defined as +, -, ^, or ~^.  For + and -, I transformed
  the output delta to a Gray code (a^(a>>1)) so a string of 1's (as
  is commonly produced by subtraction) look like a single 1-bit
  difference.
* the base values were pseudorandom, all zero but one bit set, or 
  all zero plus a counter that starts at zero.

Some k values for my "a-=c; a^=rot(c,k); c+=b;" arrangement that
satisfy this are
    4  6  8 16 19  4
    9 15  3 18 27 15
   14  9  3  7 17  3
Well, "9 15 3 18 27 15" didn't quite get 32 bits diffing
for "differ" defined as + with a one-bit base and a two-bit delta.  I
used http://burtleburtle.net/bob/hash/avalanche.html to choose 
the operations, constants, and arrangements of the variables.

This does not achieve avalanche.  There are input bits of (a,b,c)
that fail to affect some output bits of (a,b,c), especially of a.  The
most thoroughly mixed value is c, but it doesn't really even achieve
avalanche in c.

This allows some parallelism.  Read-after-writes are good at doubling
the number of bits affected, so the goal of mixing pulls in the opposite
direction as the goal of parallelism.  I did what I could.  Rotates
seem to cost as much as shifts on every machine I could lay my hands
on, and rotates are much kinder to the top and bottom bits, so I used
rotates.
-------------------------------------------------------------------------------
*/
#define mix(a,b,c) \
{ \
  a -= c;  a ^= rot(c, 4);  c += b; \
  b -= a;  b ^= rot(a, 6);  a += c; \
  c -= b;  c ^= rot(b, 8);  b += a; \
  a -= c;  a ^= rot(c,16);  c += b; \
  b -= a;  b ^= rot(a,19);  a += c; \
  c -= b;  c ^= rot(b, 4);  b += a; \
}

/*
-------------------------------------------------------------------------------
final -- final mixing of 3 32-bit values (a,b,c) into c

Pairs of (a,b,c) values differing in only a few bits will usually
produce values of c that look totally different.  This was tested for
* pairs that differed by one bit, by two bits, in any combination
  of top bits of (a,b,c), or in any combination of bottom bits of
  (a,b,c).
* "differ" is defined as +, -, ^, or ~^.  For + and -, I transformed
  the output delta to a Gray code (a^(a>>1)) so a string of 1's (as
  is commonly produced by subtraction) look like a single 1-bit
  difference.
* the base values were pseudorandom, all zero but one bit set, or 
  all zero plus a counter that starts at zero.

These constants passed:
 14 11 25 16 4 14 24
 12 14 25 16 4 14 24
and these came close:
  4  8 15 26 3 22 24
 10  8 15 26 3 22 24
 11  8 15 26 3 22 24
-------------------------------------------------------------------------------
*/
#define final(a,b,c) \
{ \
  c ^= b; c -= rot(b,14); \
  a ^= c; a -= rot(c,11); \
  b ^= a; b -= rot(a,25); \
  c ^= b; c -= rot(b,16); \
  a ^= c; a -= rot(c,4);  \
  b ^= a; b -= rot(a,14); \
  c ^= b; c -= rot(b,24); \
}

/*
--------------------------------------------------------------------
 This works on all machines.  To be useful, it requires
 -- that the key be an array of uint32_t's, and
 -- that the length be the number of uint32_t's in the key

 The function hashword() is identical to hashlittle() on little-endian
 machines, and identical to hashbig() on big-endian machines,
 except that the length has to be measured in uint32_ts rather than in
 bytes.  hashlittle() is more complicated than hashword() only because
 hashlittle() has to dance around fitting the key bytes into registers.
--------------------------------------------------------------------
*/
uint32_t hashword(
const uint32_t *k,                   /* the key, an array of uint32_t values */
size_t          length,               /* the length of the key, in uint32_ts */
uint32_t        initval)         /* the previous hash, or an arbitrary value */
{
  uint32_t a,b,c;

  /* Set up the internal state */
  a = b = c = 0xdeadbeef + (((uint32_t)length)<<2) + initval;

  /*------------------------------------------------- handle most of the key */
  while (length > 3)
  {
    a += k[0];
    b += k[1];
    c += k[2];
    mix(a,b,c);
    length -= 3;
    k += 3;
  }

  /*------------------------------------------- handle the last 3 uint32_t's */
  switch(length)                     /* all the case statements fall through */
  { 
  case 3 : c+=k[2];
  case 2 : b+=k[1];
  case 1 : a+=k[0];
    final(a,b,c);
  case 0:     /* case 0: nothing left to add */
    break;
  }
  /*------------------------------------------------------ report the result */
  return c;
}


/*
--------------------------------------------------------------------
hashword2() -- same as hashword(), but take two seeds and return two
32-bit values.  pc and pb must both be nonnull, and *pc and *pb must
both be initialized with seeds.  If you pass in (*pb)==0, the output 
(*pc) will be the same as the return value from hashword().
--------------------------------------------------------------------
*/
void hashword2 (
const uint32_t *k,                   /* the key, an array of uint32_t values */
size_t          length,               /* the length of the key, in uint32_ts */
uint32_t       *pc,                      /* IN: seed OUT: primary hash value */
uint32_t       *pb)               /* IN: more seed OUT: secondary hash value */
{
  uint32_t a,b,c;

  /* Set up the internal state */
  a = b = c = 0xdeadbeef + ((uint32_t)(length<<2)) + *pc;
  c += *pb;

  /*------------------------------------------------- handle most of the key */
  while (length > 3)
  {
    a += k[0];
    b += k[1];
    c += k[2];
    mix(a,b,c);
    length -= 3;
    k += 3;
  }

  /*------------------------------------------- handle the last 3 uint32_t's */
  switch(length)                     /* all the case statements fall through */
  { 
  case 3 : c+=k[2];
  case 2 : b+=k[1];
  case 1 : a+=k[0];
    final(a,b,c);
  case 0:     /* case 0: nothing left to add */
    break;
  }
  /*------------------------------------------------------ report the result */
  *pc=c; *pb=b;
}


/*
-------------------------------------------------------------------------------
hashlittle() -- hash a variable-length key into a 32-bit value
  k       : the key (the unaligned variable-length array of bytes)
  length  : the length of the key, counting by bytes
  initval : can be any 4-byte value
Returns a 32-bit value.  Every bit of the key affects every bit of
the return value.  Two keys differing by one or two bits will have
totally different hash values.

The best hash table sizes are powers of 2.  There is no need to do
mod a prime (mod is sooo slow!).  If you need less than 32 bits,
use a bitmask.  For example, if you need only 10 bits, do
  h = (h & hashmask(10));
In which case, the hash table should have hashsize(10) elements.

If you are hashing n strings (uint8_t **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = hashlittle( k[i], len[i], h);

By Bob Jenkins, 2006.  bob_jenkins@burtleburtle.net.  You may use this
code any way you wish, private, educational, or commercial.  It's free.

Use for hash table lookup, or anything where one collision in 2^^32 is
acceptable.  Do NOT use for cryptographic purposes.
-------------------------------------------------------------------------------
*/

uint32_t hashlittle( const void *key, size_t length, uint32_t initval)
{
  uint32_t a,b,c;                                          /* internal state */
  union { const void *ptr; size_t i; } u;     /* needed for Mac Powerbook G4 */

  /* Set up the internal state */
  a = b = c = 0xdeadbeef + ((uint32_t)length) + initval;

  u.ptr = key;
  if (HASH_LITTLE_ENDIAN && ((u.i & 0x3) == 0)) {
    const uint32_t *k = (const uint32_t *)key;         /* read 32-bit chunks */

    /*------ all but last block: aligned reads and affect 32 bits of (a,b,c) */
    while (length > 12)
    {
      a += k[0];
      b += k[1];
      c += k[2];
      mix(a,b,c);
      length -= 12;
      k += 3;
    }

    /*----------------------------- handle the last (probably partial) block */
    /* 
     * "k[2]&0xffffff" actually reads beyond the end of the string, but
     * then masks off the part it's not allowed to read.  Because the
     * string is aligned, the masked-off tail is in the same word as the
     * rest of the string.  Every machine with memory protection I've seen
     * does it on word boundaries, so is OK with this.  But VALGRIND will
     * still catch it and complain.  The masking trick does make the hash
     * noticably faster for short strings (like English words).
     */
#ifndef VALGRIND

    switch(length)
    {
    case 12: c+=k[2]; b+=k[1]; a+=k[0]; break;
    case 11: c+=k[2]&0xffffff; b+=k[1]; a+=k[0]; break;
    case 10: c+=k[2]&0xffff; b+=k[1]; a+=k[0]; break;
    case 9 : c+=k[2]&0xff; b+=k[1]; a+=k[0]; break;
    case 8 : b+=k[1]; a+=k[0]; break;
    case 7 : b+=k[1]&0xffffff; a+=k[0]; break;
    case 6 : b+=k[1]&0xffff; a+=k[0]; break;
    case 5 : b+=k[1]&0xff; a+=k[0]; break;
    case 4 : a+=k[0]; break;
    case 3 : a+=k[0]&0xffffff; break;
    case 2 : a+=k[0]&0xffff; break;
    case 1 : a+=k[0]&0xff; break;
    case 0 : return c;              /* zero length strings require no mixing */
    }

#else /* make valgrind happy */

    k8 = (const uint8_t *)k;
    switch(length)
    {
    case 12: c+=k[2]; b+=k[1]; a+=k[0]; break;
    case 11: c+=((uint32_t)k8[10])<<16;  /* fall through */
    case 10: c+=((uint32_t)k8[9])<<8;    /* fall through */
    case 9 : c+=k8[8];                   /* fall through */
    case 8 : b+=k[1]; a+=k[0]; break;
    case 7 : b+=((uint32_t)k8[6])<<16;   /* fall through */
    case 6 : b+=((uint32_t)k8[5])<<8;    /* fall through */
    case 5 : b+=k8[4];                   /* fall through */
    case 4 : a+=k[0]; break;
    case 3 : a+=((uint32_t)k8[2])<<16;   /* fall through */
    case 2 : a+=((uint32_t)k8[1])<<8;    /* fall through */
    case 1 : a+=k8[0]; break;
    case 0 : return c;
    }

#endif /* !valgrind */

  } else if (HASH_LITTLE_ENDIAN && ((u.i & 0x1) == 0)) {
    const uint16_t *k = (const uint16_t *)key;         /* read 16-bit chunks */
    const uint8_t  *k8;

    /*--------------- all but last block: aligned reads and different mixing */
    while (length > 12)
    {
      a += k[0] + (((uint32_t)k[1])<<16);
      b += k[2] + (((uint32_t)k[3])<<16);
      c += k[4] + (((uint32_t)k[5])<<16);
      mix(a,b,c);
      length -= 12;
      k += 6;
    }

    /*----------------------------- handle the last (probably partial) block */
    k8 = (const uint8_t *)k;
    switch(length)
    {
    case 12: c+=k[4]+(((uint32_t)k[5])<<16);
             b+=k[2]+(((uint32_t)k[3])<<16);
             a+=k[0]+(((uint32_t)k[1])<<16);
             break;
    case 11: c+=((uint32_t)k8[10])<<16;     /* fall through */
    case 10: c+=k[4];
             b+=k[2]+(((uint32_t)k[3])<<16);
             a+=k[0]+(((uint32_t)k[1])<<16);
             break;
    case 9 : c+=k8[8];                      /* fall through */
    case 8 : b+=k[2]+(((uint32_t)k[3])<<16);
             a+=k[0]+(((uint32_t)k[1])<<16);
             break;
    case 7 : b+=((uint32_t)k8[6])<<16;      /* fall through */
    case 6 : b+=k[2];
             a+=k[0]+(((uint32_t)k[1])<<16);
             break;
    case 5 : b+=k8[4];                      /* fall through */
    case 4 : a+=k[0]+(((uint32_t)k[1])<<16);
             break;
    case 3 : a+=((uint32_t)k8[2])<<16;      /* fall through */
    case 2 : a+=k[0];
             break;
    case 1 : a+=k8[0];
             break;
    case 0 : return c;                     /* zero length requires no mixing */
    }

  } else {                        /* need to read the key one byte at a time */
    const uint8_t *k = (const uint8_t *)key;

    /*--------------- all but the last block: affect some 32 bits of (a,b,c) */
    while (length > 12)
    {
      a += k[0];
      a += ((uint32_t)k[1])<<8;
      a += ((uint32_t)k[2])<<16;
      a += ((uint32_t)k[3])<<24;
      b += k[4];
      b += ((uint32_t)k[5])<<8;
      b += ((uint32_t)k[6])<<16;
      b += ((uint32_t)k[7])<<24;
      c += k[8];
      c += ((uint32_t)k[9])<<8;
      c += ((uint32_t)k[10])<<16;
      c += ((uint32_t)k[11])<<24;
      mix(a,b,c);
      length -= 12;
      k += 12;
    }

    /*-------------------------------- last block: affect all 32 bits of (c) */
    switch(length)                   /* all the case statements fall through */
    {
    case 12: c+=((uint32_t)k[11])<<24;
    case 11: c+=((uint32_t)k[10])<<16;
    case 10: c+=((uint32_t)k[9])<<8;
    case 9 : c+=k[8];
    case 8 : b+=((uint32_t)k[7])<<24;
    case 7 : b+=((uint32_t)k[6])<<16;
    case 6 : b+=((uint32_t)k[5])<<8;
    case 5 : b+=k[4];
    case 4 : a+=((uint32_t)k[3])<<24;
    case 3 : a+=((uint32_t)k[2])<<16;
    case 2 : a+=((uint32_t)k[1])<<8;
    case 1 : a+=k[0];
             break;
    case 0 : return c;
    }
  }

  final(a,b,c);
  return c;
}


/*
 * hashlittle2: return 2 32-bit hash values
 *
 * This is identical to hashlittle(), except it returns two 32-bit hash
 * values instead of just one.  This is good enough for hash table
 * lookup with 2^^64 buckets, or if you want a second hash if you're not
 * happy with the first, or if you want a probably-unique 64-bit ID for
 * the key.  *pc is better mixed than *pb, so use *pc first.  If you want
 * a 64-bit value do something like "*pc + (((uint64_t)*pb)<<32)".
 */
void hashlittle2( 
  const void *key,       /* the key to hash */
  size_t      length,    /* length of the key */
  uint32_t   *pc,        /* IN: primary initval, OUT: primary hash */
  uint32_t   *pb)        /* IN: secondary initval, OUT: secondary hash */
{
  uint32_t a,b,c;                                          /* internal state */
  union { const void *ptr; size_t i; } u;     /* needed for Mac Powerbook G4 */

  /* Set up the internal state */
  a = b = c = 0xdeadbeef + ((uint32_t)length) + *pc;
  c += *pb;

  u.ptr = key;
  if (HASH_LITTLE_ENDIAN && ((u.i & 0x3) == 0)) {
    const uint32_t *k = (const uint32_t *)key;         /* read 32-bit chunks */

    /*------ all but last block: aligned reads and affect 32 bits of (a,b,c) */
    while (length > 12)
    {
      a += k[0];
      b += k[1];
      c += k[2];
      mix(a,b,c);
      length -= 12;
      k += 3;
    }

    /*----------------------------- handle the last (probably partial) block */
    /* 
     * "k[2]&0xffffff" actually reads beyond the end of the string, but
     * then masks off the part it's not allowed to read.  Because the
     * string is aligned, the masked-off tail is in the same word as the
     * rest of the string.  Every machine with memory protection I've seen
     * does it on word boundaries, so is OK with this.  But VALGRIND will
     * still catch it and complain.  The masking trick does make the hash
     * noticably faster for short strings (like English words).
     */
#ifndef VALGRIND

    switch(length)
    {
    case 12: c+=k[2]; b+=k[1]; a+=k[0]; break;
    case 11: c+=k[2]&0xffffff; b+=k[1]; a+=k[0]; break;
    case 10: c+=k[2]&0xffff; b+=k[1]; a+=k[0]; break;
    case 9 : c+=k[2]&0xff; b+=k[1]; a+=k[0]; break;
    case 8 : b+=k[1]; a+=k[0]; break;
    case 7 : b+=k[1]&0xffffff; a+=k[0]; break;
    case 6 : b+=k[1]&0xffff; a+=k[0]; break;
    case 5 : b+=k[1]&0xff; a+=k[0]; break;
    case 4 : a+=k[0]; break;
    case 3 : a+=k[0]&0xffffff; break;
    case 2 : a+=k[0]&0xffff; break;
    case 1 : a+=k[0]&0xff; break;
    case 0 : *pc=c; *pb=b; return;  /* zero length strings require no mixing */
    }

#else /* make valgrind happy */

    k8 = (const uint8_t *)k;
    switch(length)
    {
    case 12: c+=k[2]; b+=k[1]; a+=k[0]; break;
    case 11: c+=((uint32_t)k8[10])<<16;  /* fall through */
    case 10: c+=((uint32_t)k8[9])<<8;    /* fall through */
    case 9 : c+=k8[8];                   /* fall through */
    case 8 : b+=k[1]; a+=k[0]; break;
    case 7 : b+=((uint32_t)k8[6])<<16;   /* fall through */
    case 6 : b+=((uint32_t)k8[5])<<8;    /* fall through */
    case 5 : b+=k8[4];                   /* fall through */
    case 4 : a+=k[0]; break;
    case 3 : a+=((uint32_t)k8[2])<<16;   /* fall through */
    case 2 : a+=((uint32_t)k8[1])<<8;    /* fall through */
    case 1 : a+=k8[0]; break;
    case 0 : *pc=c; *pb=b; return;  /* zero length strings require no mixing */
    }

#endif /* !valgrind */

  } else if (HASH_LITTLE_ENDIAN && ((u.i & 0x1) == 0)) {
    const uint16_t *k = (const uint16_t *)key;         /* read 16-bit chunks */
    const uint8_t  *k8;

    /*--------------- all but last block: aligned reads and different mixing */
    while (length > 12)
    {
      a += k[0] + (((uint32_t)k[1])<<16);
      b += k[2] + (((uint32_t)k[3])<<16);
      c += k[4] + (((uint32_t)k[5])<<16);
      mix(a,b,c);
      length -= 12;
      k += 6;
    }

    /*----------------------------- handle the last (probably partial) block */
    k8 = (const uint8_t *)k;
    switch(length)
    {
    case 12: c+=k[4]+(((uint32_t)k[5])<<16);
             b+=k[2]+(((uint32_t)k[3])<<16);
             a+=k[0]+(((uint32_t)k[1])<<16);
             break;
    case 11: c+=((uint32_t)k8[10])<<16;     /* fall through */
    case 10: c+=k[4];
             b+=k[2]+(((uint32_t)k[3])<<16);
             a+=k[0]+(((uint32_t)k[1])<<16);
             break;
    case 9 : c+=k8[8];                      /* fall through */
    case 8 : b+=k[2]+(((uint32_t)k[3])<<16);
             a+=k[0]+(((uint32_t)k[1])<<16);
             break;
    case 7 : b+=((uint32_t)k8[6])<<16;      /* fall through */
    case 6 : b+=k[2];
             a+=k[0]+(((uint32_t)k[1])<<16);
             break;
    case 5 : b+=k8[4];                      /* fall through */
    case 4 : a+=k[0]+(((uint32_t)k[1])<<16);
             break;
    case 3 : a+=((uint32_t)k8[2])<<16;      /* fall through */
    case 2 : a+=k[0];
             break;
    case 1 : a+=k8[0];
             break;
    case 0 : *pc=c; *pb=b; return;  /* zero length strings require no mixing */
    }

  } else {                        /* need to read the key one byte at a time */
    const uint8_t *k = (const uint8_t *)key;

    /*--------------- all but the last block: affect some 32 bits of (a,b,c) */
    while (length > 12)
    {
      a += k[0];
      a += ((uint32_t)k[1])<<8;
      a += ((uint32_t)k[2])<<16;
      a += ((uint32_t)k[3])<<24;
      b += k[4];
      b += ((uint32_t)k[5])<<8;
      b += ((uint32_t)k[6])<<16;
      b += ((uint32_t)k[7])<<24;
      c += k[8];
      c += ((uint32_t)k[9])<<8;
      c += ((uint32_t)k[10])<<16;
      c += ((uint32_t)k[11])<<24;
      mix(a,b,c);
      length -= 12;
      k += 12;
    }

    /*-------------------------------- last block: affect all 32 bits of (c) */
    switch(length)                   /* all the case statements fall through */
    {
    case 12: c+=((uint32_t)k[11])<<24;
    case 11: c+=((uint32_t)k[10])<<16;
    case 10: c+=((uint32_t)k[9])<<8;
    case 9 : c+=k[8];
    case 8 : b+=((uint32_t)k[7])<<24;
    case 7 : b+=((uint32_t)k[6])<<16;
    case 6 : b+=((uint32_t)k[5])<<8;
    case 5 : b+=k[4];
    case 4 : a+=((uint32_t)k[3])<<24;
    case 3 : a+=((uint32_t)k[2])<<16;
    case 2 : a+=((uint32_t)k[1])<<8;
    case 1 : a+=k[0];
             break;
    case 0 : *pc=c; *pb=b; return;  /* zero length strings require no mixing */
    }
  }

  final(a,b,c);
  *pc=c; *pb=b;
}



/*
 * hashbig():
 * This is the same as hashword() on big-endian machines.  It is different
 * from hashlittle() on all machines.  hashbig() takes advantage of
 * big-endian byte ordering. 
 */
uint32_t hashbig( const void *key, size_t length, uint32_t initval)
{
  uint32_t a,b,c;
  union { const void *ptr; size_t i; } u; /* to cast key to (size_t) happily */

  /* Set up the internal state */
  a = b = c = 0xdeadbeef + ((uint32_t)length) + initval;

  u.ptr = key;
  if (HASH_BIG_ENDIAN && ((u.i & 0x3) == 0)) {
    const uint32_t *k = (const uint32_t *)key;         /* read 32-bit chunks */

    /*------ all but last block: aligned reads and affect 32 bits of (a,b,c) */
    while (length > 12)
    {
      a += k[0];
      b += k[1];
      c += k[2];
      mix(a,b,c);
      length -= 12;
      k += 3;
    }

    /*----------------------------- handle the last (probably partial) block */
    /* 
     * "k[2]<<8" actually reads beyond the end of the string, but
     * then shifts out the part it's not allowed to read.  Because the
     * string is aligned, the illegal read is in the same word as the
     * rest of the string.  Every machine with memory protection I've seen
     * does it on word boundaries, so is OK with this.  But VALGRIND will
     * still catch it and complain.  The masking trick does make the hash
     * noticably faster for short strings (like English words).
     */
#ifndef VALGRIND

    switch(length)
    {
    case 12: c+=k[2]; b+=k[1]; a+=k[0]; break;
    case 11: c+=k[2]&0xffffff00; b+=k[1]; a+=k[0]; break;
    case 10: c+=k[2]&0xffff0000; b+=k[1]; a+=k[0]; break;
    case 9 : c+=k[2]&0xff000000; b+=k[1]; a+=k[0]; break;
    case 8 : b+=k[1]; a+=k[0]; break;
    case 7 : b+=k[1]&0xffffff00; a+=k[0]; break;
    case 6 : b+=k[1]&0xffff0000; a+=k[0]; break;
    case 5 : b+=k[1]&0xff000000; a+=k[0]; break;
    case 4 : a+=k[0]; break;
    case 3 : a+=k[0]&0xffffff00; break;
    case 2 : a+=k[0]&0xffff0000; break;
    case 1 : a+=k[0]&0xff000000; break;
    case 0 : return c;              /* zero length strings require no mixing */
    }

#else  /* make valgrind happy */

    k8 = (const uint8_t *)k;
    switch(length)                   /* all the case statements fall through */
    {
    case 12: c+=k[2]; b+=k[1]; a+=k[0]; break;
    case 11: c+=((uint32_t)k8[10])<<8;  /* fall through */
    case 10: c+=((uint32_t)k8[9])<<16;  /* fall through */
    case 9 : c+=((uint32_t)k8[8])<<24;  /* fall through */
    case 8 : b+=k[1]; a+=k[0]; break;
    case 7 : b+=((uint32_t)k8[6])<<8;   /* fall through */
    case 6 : b+=((uint32_t)k8[5])<<16;  /* fall through */
    case 5 : b+=((uint32_t)k8[4])<<24;  /* fall through */
    case 4 : a+=k[0]; break;
    case 3 : a+=((uint32_t)k8[2])<<8;   /* fall through */
    case 2 : a+=((uint32_t)k8[1])<<16;  /* fall through */
    case 1 : a+=((uint32_t)k8[0])<<24; break;
    case 0 : return c;
    }

#endif /* !VALGRIND */

  } else {                        /* need to read the key one byte at a time */
    const uint8_t *k = (const uint8_t *)key;

    /*--------------- all but the last block: affect some 32 bits of (a,b,c) */
    while (length > 12)
    {
      a += ((uint32_t)k[0])<<24;
      a += ((uint32_t)k[1])<<16;
      a += ((uint32_t)k[2])<<8;
      a += ((uint32_t)k[3]);
      b += ((uint32_t)k[4])<<24;
      b += ((uint32_t)k[5])<<16;
      b += ((uint32_t)k[6])<<8;
      b += ((uint32_t)k[7]);
      c += ((uint32_t)k[8])<<24;
      c += ((uint32_t)k[9])<<16;
      c += ((uint32_t)k[10])<<8;
      c += ((uint32_t)k[11]);
      mix(a,b,c);
      length -= 12;
      k += 12;
    }

    /*-------------------------------- last block: affect all 32 bits of (c) */
    switch(length)                   /* all the case statements fall through */
    {
    case 12: c+=k[11];
    case 11: c+=((uint32_t)k[10])<<8;
    case 10: c+=((uint32_t)k[9])<<16;
    case 9 : c+=((uint32_t)k[8])<<24;
    case 8 : b+=k[7];
    case 7 : b+=((uint32_t)k[6])<<8;
    case 6 : b+=((uint32_t)k[5])<<16;
    case 5 : b+=((uint32_t)k[4])<<24;
    case 4 : a+=k[3];
    case 3 : a+=((uint32_t)k[2])<<8;
    case 2 : a+=((uint32_t)k[1])<<16;
    case 1 : a+=((uint32_t)k[0])<<24;
             break;
    case 0 : return c;
    }
  }

  final(a,b,c);
  return c;
}




