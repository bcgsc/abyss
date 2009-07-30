#include "PackedSeq.h"
#include "HashFunction.h"
#include "Options.h"
#include "Sequence.h"
#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace std;

//
// Default constructor
//
PackedSeq::PackedSeq() : m_length(0), m_flags(0)
{
	m_multiplicity[SENSE] = 1;
	m_multiplicity[ANTISENSE] = 0;
}

//
// Construct a sequence from a String-based sequence
//
PackedSeq::PackedSeq(const Sequence& seq)
	: m_length(seq.length()), m_flags(0)
{
	memset(m_seq, 0, NUM_BYTES);
	assert(m_length <= MAX_KMER);
	const char* p = seq.data();
	for(unsigned i = 0; i < m_length; i++)
		setBaseCode(m_seq, i, baseToCode(*p++));
	m_multiplicity[SENSE] = 1;
	m_multiplicity[ANTISENSE] = 0;
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
	
	memcpy(buffer + offset, m_multiplicity, sizeof(m_multiplicity));
	offset += sizeof(m_multiplicity);	
	
	memcpy(buffer + offset, &m_extRecord, sizeof(m_extRecord));
	offset += sizeof(m_extRecord);
	
	assert(offset == serialSize());

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
	
	memcpy(m_multiplicity, buffer + offset, sizeof(m_multiplicity));
	offset += sizeof(m_multiplicity);	
	
	memcpy(&m_extRecord, buffer + offset, sizeof(m_extRecord));
	offset += sizeof(m_extRecord);
	
	assert(offset == serialSize());
	
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
	unsigned numBytes = getNumCodingBytes(m_length);

	// copy the sequence over
	memcpy(m_seq, other.m_seq, numBytes);
	
	m_flags = other.m_flags;
	
	
	m_extRecord.dir[SENSE] = other.m_extRecord.dir[SENSE];
	m_extRecord.dir[ANTISENSE] = other.m_extRecord.dir[ANTISENSE];
	m_multiplicity[SENSE] = other.m_multiplicity[SENSE];
	m_multiplicity[ANTISENSE] = other.m_multiplicity[ANTISENSE];

	return *this;
}

//
// Compare two sequences.
// 
int PackedSeq::compare(const PackedSeq& other) const
{
	if (m_length == 0 || other.m_length == 0)
		return (int)m_length - other.m_length;
	assert(m_length == other.m_length);
	unsigned nbytes = getNumCodingBytes(m_length);
	return memcmp(m_seq, other.m_seq, nbytes);
}

//
// Get the length of the sequence
//
unsigned PackedSeq::getSequenceLength() const
{
	return m_length;
}

//
// Get the number of coding bytes
// 
unsigned PackedSeq::getNumCodingBytes(unsigned seqLength)
{
	return (seqLength + 3) / 4;
}

/** Compute a hash-like value of the packed sequence over the first 16
 * bases and the reverse complement of the last 16 bases
 * The reverse complement of the last 16 bases is used so that a
 * sequence and its reverse complement will hash to the same value.
 * @todo make this faster
 */
unsigned PackedSeq::getCode() const
{
	/* At k=19, this hash function always returns a positive number
	 * due to the sequence and its reverse complement overlapping when
	 * the xor is calculated. A more general solution is needed. */
	const unsigned NUM_BYTES = m_length < 20 ? m_length/8 : 4;
	PackedSeq rc = *this;
	rc.reverseComplement();

	const unsigned prime = 101;
	unsigned sum = 0;
	for (unsigned i = 0; i < NUM_BYTES; i++)
		sum = prime * sum + (m_seq[i] ^ rc.m_seq[i]);
	return sum;
}

size_t PackedSeq::getHashCode() const
{
	// Hash on the numbytes - 1. This is to avoid getting different hash values for the same sequence for n % 4 != 0 sequences
	int code = hashlittle(m_seq, getNumCodingBytes(m_length) - 1, 131);

	return code;
}

/** Return the string representation of this sequence. */
Sequence PackedSeq::decode() const
{
	Sequence s;
	s.reserve(m_length);
	for (unsigned i = 0; i < m_length; i++)
		s.push_back(codeToBase(getBaseCode(i)));
	return s;
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

#define SEQ_WORDS ((PackedSeq::NUM_BYTES + 7)/8)
#define SEQ_BITS (64 * SEQ_WORDS)
#define SEQ_FULL_WORDS ((int)PackedSeq::NUM_BYTES/8)
#define SEQ_ODD_BYTES (PackedSeq::NUM_BYTES - 8*SEQ_FULL_WORDS)

struct Seq {
	uint64_t x[SEQ_WORDS];
};

/** Load with appropriate endianness for shifting. */
static Seq load(const uint8_t *src)
{
	Seq seq;
	uint64_t *px = seq.x;
	const uint8_t *p = src;
	for (int i = 0; i < SEQ_FULL_WORDS; i++) {
		*px++ = (uint64_t)p[0] << 56
			| (uint64_t)p[1] << 48
			| (uint64_t)p[2] << 40
			| (uint64_t)p[3] << 32
			| (uint64_t)p[4] << 24
			| (uint64_t)p[5] << 16
			| (uint64_t)p[6] << 8
			| (uint64_t)p[7] << 0;
		p += 8;
	}
	if (SEQ_ODD_BYTES > 0) {
		uint64_t x = 0;
		if (SEQ_ODD_BYTES > 0) x |= (uint64_t)p[0] << 56;
		if (SEQ_ODD_BYTES > 1) x |= (uint64_t)p[1] << 48;
		if (SEQ_ODD_BYTES > 2) x |= (uint64_t)p[2] << 40;
		if (SEQ_ODD_BYTES > 3) x |= (uint64_t)p[3] << 32;
		if (SEQ_ODD_BYTES > 4) x |= (uint64_t)p[4] << 24;
		if (SEQ_ODD_BYTES > 5) x |= (uint64_t)p[5] << 16;
		if (SEQ_ODD_BYTES > 6) x |= (uint64_t)p[6] << 8;
		if (SEQ_ODD_BYTES > 7) x |= (uint64_t)p[7] << 0;
		*px = x;
	}
	return seq;
}

static void store(uint8_t *dest, Seq seq)
{
	const uint64_t *px = seq.x;
	for (int i = 0; i < SEQ_FULL_WORDS; i++) {
		dest[0] = *px >> 56;
		dest[1] = *px >> 48;
		dest[2] = *px >> 40;
		dest[3] = *px >> 32;
		dest[4] = *px >> 24;
		dest[5] = *px >> 16;
		dest[6] = *px >> 8;
		dest[7] = *px >> 0;
		dest += 8;
		px++;
	}
	if (SEQ_ODD_BYTES > 0) dest[0] = *px >> 56;
	if (SEQ_ODD_BYTES > 1) dest[1] = *px >> 48;
	if (SEQ_ODD_BYTES > 2) dest[2] = *px >> 40;
	if (SEQ_ODD_BYTES > 3) dest[3] = *px >> 32;
	if (SEQ_ODD_BYTES > 4) dest[4] = *px >> 24;
	if (SEQ_ODD_BYTES > 5) dest[5] = *px >> 16;
	if (SEQ_ODD_BYTES > 6) dest[6] = *px >> 8;
	if (SEQ_ODD_BYTES > 7) dest[7] = *px >> 0;
}

/**
 * Reverse the bytes by storing them in the reverse order of
 * loading, and reverse the words in the same fasion.
 */
static void storeReverse(uint8_t *dest, Seq seq)
{
	uint64_t *d = (uint64_t*)dest;
	uint64_t *px = &seq.x[SEQ_WORDS-1];
	for (int i = 0; i < SEQ_FULL_WORDS; i++)
		*d++ = *px--;
	dest = (uint8_t*)d;
	if (SEQ_ODD_BYTES > 0) dest[0] = *px >> 0;
	if (SEQ_ODD_BYTES > 1) dest[1] = *px >> 8;
	if (SEQ_ODD_BYTES > 2) dest[2] = *px >> 16;
	if (SEQ_ODD_BYTES > 3) dest[3] = *px >> 24;
	if (SEQ_ODD_BYTES > 4) dest[4] = *px >> 32;
	if (SEQ_ODD_BYTES > 5) dest[5] = *px >> 40;
	if (SEQ_ODD_BYTES > 6) dest[6] = *px >> 48;
	if (SEQ_ODD_BYTES > 7) dest[7] = *px >> 56;
}

/** Shift right by the specified number of bits. */
static void shiftRight(Seq *pseq, uint8_t n)
{
	if (n == 0)
		return;
#if MAX_KMER <= 32
	pseq->x[0] >>= n;
#elif MAX_KMER <= 64
	uint64_t x0 = pseq->x[0], x1 = pseq->x[1];
	if (n < 64) {
		pseq->x[0] = x0 >> n;
		pseq->x[1] = x1 >> n | x0 << (64 - n);
	} else {
		pseq->x[0] = 0;
		pseq->x[1] = x0 >> (n - 64);
	}
#elif MAX_KMER <= 96
	uint64_t x0 = pseq->x[0], x1 = pseq->x[1], x2 = pseq->x[2];
	if (n < 64) {
		pseq->x[0] = x0 >> n;
		pseq->x[1] = x1 >> n | x0 << (64 - n);
		pseq->x[2] = x2 >> n | x1 << (64 - n);
	} else if (n == 64) {
		pseq->x[0] = 0;
		pseq->x[1] = x0;
		pseq->x[2] = x1;
	} else if (n < 128) {
		n -= 64;
		pseq->x[0] = 0;
		pseq->x[1] = x0 >> n;
		pseq->x[2] = x1 >> n | x0 << (64 - n);
	} else {
		n -= 128;
		pseq->x[0] = 0;
		pseq->x[1] = 0;
		pseq->x[2] = x0 >> n;
	}
#endif
}

/** Shift left by the specified number of bits. */
static void shiftLeft(Seq *pseq, uint8_t n)
{
	if (n == 0)
		return;
#if MAX_KMER <= 32
	pseq->x[0] <<= n;
#elif MAX_KMER <= 64
	uint64_t x0 = pseq->x[0], x1 = pseq->x[1];
	if (n < 64) {
		pseq->x[0] = x0 << n | x1 >> (64 - n);
		pseq->x[1] = x1 << n;
	} else {
		pseq->x[0] = x1 << (n - 64);
		pseq->x[1] = 0;
	}
#elif MAX_KMER <= 96
	uint64_t x0 = pseq->x[0], x1 = pseq->x[1], x2 = pseq->x[2];
	if (n < 64) {
		pseq->x[0] = x0 << n | x1 >> (64 - n);
		pseq->x[1] = x1 << n | x2 >> (64 - n);
		pseq->x[2] = x2 << n;
	} else if (n == 64) {
		pseq->x[0] = x1;
		pseq->x[1] = x2;
		pseq->x[2] = 0;
	} else if (n < 128) {
		n -= 64;
		pseq->x[0] = x1 << n | x2 >> (64 - n);
		pseq->x[1] = x2 << n;
		pseq->x[2] = 0;
	} else {
		n -= 128;
		pseq->x[0] = x2 << n;
		pseq->x[1] = 0;
		pseq->x[2] = 0;
	}
#endif
}

//
// Return a subsequence of this sequence.
//
PackedSeq PackedSeq::subseq(unsigned start, unsigned len) const
{
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
	shiftRight(&seq, SEQ_BITS - 2*(len + start));
	shiftLeft(&seq, SEQ_BITS - 2*len);
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
	if (!opt::colourSpace)
		for (unsigned i = 0; i < SEQ_WORDS; i++)
			seq.x[i] = ~seq.x[i];

	// Shift the bits flush to the right of the double word.
	shiftRight(&seq, SEQ_BITS - 2*m_length);

	storeReverse((uint8_t*)m_seq, seq);

	// Reverse the pairs of bits withing a byte.
	unsigned numBytes = getNumCodingBytes(m_length);
	for (unsigned i = 0; i < numBytes; i++)
		m_seq[i] = swapBases[(uint8_t)m_seq[i]];
}

uint8_t PackedSeq::shift(extDirection dir, uint8_t base)
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

void PackedSeq::setLastBase(extDirection dir, uint8_t base)
{
	setBaseCode(m_seq, dir == SENSE ? m_length - 1 : 0, base);
}

uint8_t PackedSeq::shiftAppend(uint8_t base)
{
	// shift the sequence left and append a new base to the end
	unsigned numBytes = getNumCodingBytes(m_length);
	uint8_t shiftIn = base;
	// starting from the last byte, shift the new base in and get the captured base
	for(int i = numBytes - 1; i >= 0; i--)
	{
		// calculate the index
		// if this is the last byte, use 
		unsigned index = (unsigned)i == numBytes - 1
			? seqIndexToBaseIndex(m_length - 1) : 3;
		shiftIn = leftShiftByte(m_seq, i, index, shiftIn);
	}
	
	// return the base shifted out of the first byte
	return shiftIn;
}

uint8_t PackedSeq::shiftPrepend(uint8_t base)
{
	// shift the sequence right and append a new base to the end
	unsigned numBytes = getNumCodingBytes(m_length);

	unsigned lastBaseByte = seqIndexToByteNumber(m_length - 1);
	unsigned lastBaseIndex = seqIndexToBaseIndex(m_length - 1);

	// save the last base (which gets shifted out)
	uint8_t lastBase = getBaseCode(m_seq,
			lastBaseByte, lastBaseIndex);
	// Zero the last base, which is required by compare.
	setBaseCode(m_seq, lastBaseByte, lastBaseIndex, 0);

	uint8_t shiftIn = base;
	// starting from the last byte, shift the new base in and get the captured base
	for(unsigned i = 0; i <= numBytes - 1; i++)
	{
		// index is always zero
		unsigned index = 0;
		shiftIn = rightShiftByte(m_seq, i, index, shiftIn);
	}
		
	return lastBase;	
}

uint8_t PackedSeq::leftShiftByte(char* pSeq,
		unsigned byteNum, unsigned index, uint8_t base)
{
	// save the first base
	uint8_t outBase = (pSeq[byteNum] >> 6) & 0x3;
	
	// shift left one position
	pSeq[byteNum] <<= 2;
	
	// Set the new base
	setBaseCode(pSeq, byteNum, index, base);

	return outBase;
}

uint8_t PackedSeq::rightShiftByte(char* pSeq,
		unsigned byteNum, unsigned index, uint8_t base)
{
	// save the last base
	uint8_t outBase = pSeq[byteNum] & 0x3;
	
	// shift right one position
	pSeq[byteNum] >>= 2;
	
	// add the new base
	setBaseCode(pSeq, byteNum, index, base);
	
	return outBase;
}

void PackedSeq::setBaseExtension(extDirection dir, uint8_t base)
{
	m_extRecord.dir[dir].setBase(base);
}

void PackedSeq::clearAllExtensions(extDirection dir)
{
	m_extRecord.dir[dir].clear();
}

void PackedSeq::clearExtension(extDirection dir, uint8_t base)
{
	m_extRecord.dir[dir].clearBase(base);
}

bool PackedSeq::hasExtension(extDirection dir) const
{
	return m_extRecord.dir[dir].hasExtension();
}

bool PackedSeq::isAmbiguous(extDirection dir) const
{
	return m_extRecord.dir[dir].isAmbiguous();
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
void PackedSeq::setBaseCode(char* pSeq,
		unsigned seqIndex, uint8_t base)
{
	return setBaseCode(pSeq,
			seqIndexToByteNumber(seqIndex),
			seqIndexToBaseIndex(seqIndex),
			base);
}

//
//Set a base by byte number/ sub index
// beware, this does not check for out of bounds access
//
void PackedSeq::setBaseCode(char* pSeq,
		unsigned byteNum, unsigned index, uint8_t base)
{
	// shift the value into position
	unsigned shiftValue = 2*(3 - index);
	base <<= shiftValue;

	// clear the value
	uint8_t mask = 0x3;
	mask <<= shiftValue;
	mask = ~mask;
	pSeq[byteNum] &= mask;

	// set the appropriate value with an OR
	pSeq[byteNum] |= base;
}

//
// get a base code by the index [0, length)
//
uint8_t PackedSeq::getBaseCode(unsigned seqIndex) const
{
	assert(seqIndex < m_length);
	return getBaseCode(m_seq,
			seqIndexToByteNumber(seqIndex),
			seqIndexToBaseIndex(seqIndex));
}

uint8_t PackedSeq::getLastBaseChar() const
{
	return codeToBase(getBaseCode(m_length - 1));
}

//
// get a base code by the byte number and sub index
//
uint8_t PackedSeq::getBaseCode(const char* pSeq,
		unsigned byteNum, unsigned index)
{
	unsigned shiftLen = 2 * (3 - index);
	return (pSeq[byteNum] >> shiftLen) & 0x3;
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

PackedSeq reverseComplement(const PackedSeq& seq)
{
	PackedSeq rc(seq);
	rc.reverseComplement();
	return rc;	
}

/** Return true if this sequence is a palindrome. */
bool PackedSeq::isPalindrome() const
{
	return m_length % 2 == 1 ? false
		: *this == ::reverseComplement(*this);
}

/** Return true if the length k-1 subsequence is a palindrome. */
bool PackedSeq::isPalindrome(extDirection dir) const
{
	if (m_length % 2 == 0)
		return false;
	PackedSeq seq(*this);
	if (dir == SENSE)
		seq.shiftAppend(0);
	else
		seq.setLastBase(SENSE, 0);
	seq.m_length--;
	return seq.isPalindrome();
}
