#include "config.h"
#include "Kmer.h"
#include "Common/Options.h"
#include "HashFunction.h"
#include <algorithm>
#include <cstring>

using namespace std;

/** The size of a k-mer. This variable is static and is shared by all
 * instances. */
unsigned Kmer::s_length;

/** The size of a k-mer in bytes. */
unsigned Kmer::s_bytes;

static unsigned seqIndexToByteNumber(unsigned seqIndex);
static unsigned seqIndexToBaseIndex(unsigned seqIndex);
static uint8_t getBaseCode(const char* pSeq,
		unsigned byteNum, unsigned index);
static void setBaseCode(char* pSeq,
		unsigned byteNum, unsigned index, uint8_t base);

/** Construct a k-mer from a string. */
Kmer::Kmer(const Sequence& seq)
{
	assert(seq.length() == s_length);
	memset(m_seq, 0, NUM_BYTES);
	const char* p = seq.data();
	for (unsigned i = 0; i < s_length; i++)
		set(i, baseToCode(*p++));
}

/** Compare two k-mer. */
int Kmer::compare(const Kmer& other) const
{
	return memcmp(m_seq, other.m_seq, bytes());
}

/** Compute a hash-like value of the packed sequence over the first 16
 * bases and the reverse complement of the last 16 bases
 * The reverse complement of the last 16 bases is used so that a
 * sequence and its reverse complement will hash to the same value.
 * @todo make this faster
 */
unsigned Kmer::getCode() const
{
	/* At k=19, this hash function always returns an even number due
	 * to the sequence and its reverse complement overlapping when the
	 * xor is calculated. A more general solution is needed. */
	const unsigned NUM_BYTES
		= s_length < 8 ? 1
		: s_length < 20 ? s_length/8
		: 4;
	Kmer rc = *this;
	rc.reverseComplement();

	const unsigned prime = 101;
	unsigned sum = 0;
	for (unsigned i = 0; i < NUM_BYTES; i++)
		sum = prime * sum + (m_seq[i] ^ rc.m_seq[i]);
	return sum;
}

size_t Kmer::getHashCode() const
{
	// Hash numbytes - 1 to avoid getting different hash values for
	// the same sequence for n % 4 != 0 sequences.
	return hashmem(m_seq, bytes() - 1);
}

/** Return the string representation of this sequence. */
Sequence Kmer::str() const
{
	Sequence s;
	s.reserve(s_length);
	for (unsigned i = 0; i < s_length; i++)
		s.push_back(codeToBase(at(i)));
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

#if MAX_KMER > 96
# include <bitset>
typedef bitset<2 * MAX_KMER> Seq;
#else
# define SEQ_WORDS ((Kmer::NUM_BYTES + 7)/8)
# define SEQ_BITS (64 * SEQ_WORDS)
# define SEQ_FULL_WORDS ((int)Kmer::NUM_BYTES/8)
# define SEQ_ODD_BYTES (Kmer::NUM_BYTES - 8*SEQ_FULL_WORDS)

/** A sequence of bits of length SEQ_BITS. */
struct Seq {
	uint64_t x[SEQ_WORDS];

	/** Return the number of bits in this sequence. */
	static unsigned size() { return SEQ_BITS; }

	/** Flip all the bits. */
	void flip()
	{
		for (unsigned i = 0; i < SEQ_WORDS; i++)
			x[i] = ~x[i];
	}

	/** Shift right by the specified number of bits. */
	void operator>>=(uint8_t n)
	{
		if (n == 0)
			return;
#if MAX_KMER <= 32
		x[0] >>= n;
#elif MAX_KMER <= 64
		uint64_t x0 = x[0], x1 = x[1];
		if (n < 64) {
			x[0] = x0 >> n;
			x[1] = x1 >> n | x0 << (64 - n);
		} else {
			x[0] = 0;
			x[1] = x0 >> (n - 64);
		}
#elif MAX_KMER <= 96
		uint64_t x0 = x[0], x1 = x[1], x2 = x[2];
		if (n < 64) {
			x[0] = x0 >> n;
			x[1] = x1 >> n | x0 << (64 - n);
			x[2] = x2 >> n | x1 << (64 - n);
		} else if (n == 64) {
			x[0] = 0;
			x[1] = x0;
			x[2] = x1;
		} else if (n < 128) {
			n -= 64;
			x[0] = 0;
			x[1] = x0 >> n;
			x[2] = x1 >> n | x0 << (64 - n);
		} else {
			n -= 128;
			x[0] = 0;
			x[1] = 0;
			x[2] = x0 >> n;
		}
#else
# error
#endif
	}
};
#endif

/** Load with appropriate endianness for shifting. */
static Seq load(const uint8_t *src)
{
	Seq seq;
#if MAX_KMER > 96
# if WORDS_BIGENDIAN
	const size_t *s = reinterpret_cast<const size_t*>(src);
	size_t *d = reinterpret_cast<size_t*>(&seq + 1);
	copy(s, s + Kmer::NUM_BYTES/sizeof(size_t), reverse_iterator<size_t*>(d));
# else
	uint8_t *d = reinterpret_cast<uint8_t*>(&seq);
	memcpy(d, src, sizeof seq);
	reverse(d, d + sizeof seq);
# endif
#else
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
#endif
	return seq;
}

/**
 * Reverse the bytes by storing them in the reverse order of
 * loading, and reverse the words in the same fashion.
 */
static void storeReverse(uint8_t *dest, const Seq seq)
{
#if MAX_KMER > 96
# if WORDS_BIGENDIAN
	const size_t *s = reinterpret_cast<const size_t*>(&seq);
	size_t *d = reinterpret_cast<size_t*>(dest);
	copy(s, s + Kmer::NUM_BYTES/sizeof(size_t),
			reverse_iterator<size_t*>(d +  Kmer::NUM_BYTES/sizeof(size_t)));
	reverse(dest, dest + Kmer::NUM_BYTES);
# else
	memcpy(dest, &seq, Kmer::NUM_BYTES);
# endif
#else
	const uint64_t *px = &seq.x[SEQ_WORDS-1];
# if WORDS_BIGENDIAN
	for (int i = 0; i < SEQ_FULL_WORDS; i++) {
		dest[0] = *px >> 0;
		dest[1] = *px >> 8;
		dest[2] = *px >> 16;
		dest[3] = *px >> 24;
		dest[4] = *px >> 32;
		dest[5] = *px >> 40;
		dest[6] = *px >> 48;
		dest[7] = *px >> 56;
		dest += 8;
		px--;
	}
# else
	uint64_t *d = (uint64_t*)dest;
	for (int i = 0; i < SEQ_FULL_WORDS; i++)
		*d++ = *px--;
	dest = (uint8_t*)d;
# endif
	if (SEQ_ODD_BYTES > 0) dest[0] = *px >> 0;
	if (SEQ_ODD_BYTES > 1) dest[1] = *px >> 8;
	if (SEQ_ODD_BYTES > 2) dest[2] = *px >> 16;
	if (SEQ_ODD_BYTES > 3) dest[3] = *px >> 24;
	if (SEQ_ODD_BYTES > 4) dest[4] = *px >> 32;
	if (SEQ_ODD_BYTES > 5) dest[5] = *px >> 40;
	if (SEQ_ODD_BYTES > 6) dest[6] = *px >> 48;
	if (SEQ_ODD_BYTES > 7) dest[7] = *px >> 56;
#endif
}

/** Reverse-complement this sequence. */
void Kmer::reverseComplement()
{
	Seq seq = load((uint8_t*)m_seq);

	// Complement the bits.
	if (!opt::colourSpace)
		seq.flip();

	// Shift the bits flush to the right of the double word.
	seq >>= seq.size() - 2*s_length;

	storeReverse((uint8_t*)m_seq, seq);

	// Reverse the pairs of bits withing a byte.
	unsigned numBytes = bytes();
	for (unsigned i = 0; i < numBytes; i++)
		m_seq[i] = swapBases[(uint8_t)m_seq[i]];
}

bool Kmer::isCanonical() const
{
	for (unsigned i = 0, j = s_length - 1;
		i < s_length / 2 + s_length % 2; i++, j--) {
		uint8_t base = getBaseCode(m_seq,
			seqIndexToByteNumber(i), seqIndexToBaseIndex(i));
		uint8_t rcBase = 0x3 & ~getBaseCode(m_seq,
			seqIndexToByteNumber(j), seqIndexToBaseIndex(j));
		if (base == rcBase)
			continue;
		return rcBase > base;
	}
	return true;
}

void Kmer::canonicalize()
{
	if (!isCanonical())
		reverseComplement();
}

void Kmer::setLastBase(extDirection dir, uint8_t base)
{
	set(dir == SENSE ? s_length - 1 : 0, base);
}

/** Shift the sequence left and append a new base to the end.
 * @return the base shifted out
 */
uint8_t Kmer::shiftAppend(uint8_t base)
{
	unsigned numBytes = bytes();
	uint8_t shiftIn = base;
	for(int i = numBytes - 1; i >= 0; i--)
	{
		unsigned index = (unsigned)i == numBytes - 1
			? seqIndexToBaseIndex(s_length - 1) : 3;
		shiftIn = leftShiftByte(m_seq, i, index, shiftIn);
	}
	return shiftIn;
}

/** Shift the sequence right and prepend a new base at the front.
 * @return the base shifted out
 */
uint8_t Kmer::shiftPrepend(uint8_t base)
{
	unsigned numBytes = bytes();

	unsigned lastBaseByte = seqIndexToByteNumber(s_length - 1);
	unsigned lastBaseIndex = seqIndexToBaseIndex(s_length - 1);

	// save the last base (which gets shifted out)
	uint8_t lastBase = getBaseCode(m_seq,
			lastBaseByte, lastBaseIndex);
	// Zero the last base, which is required by compare.
	setBaseCode(m_seq, lastBaseByte, lastBaseIndex, 0);

	uint8_t shiftIn = base;
	for(unsigned i = 0; i <= numBytes - 1; i++)
	{
		// index is always zero
		unsigned index = 0;
		shiftIn = rightShiftByte(m_seq, i, index, shiftIn);
	}
	return lastBase;
}

uint8_t Kmer::leftShiftByte(char* pSeq,
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

uint8_t Kmer::rightShiftByte(char* pSeq,
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

//Set a base by byte number/ sub index
// beware, this does not check for out of bounds access
static void setBaseCode(char* pSeq,
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

/** Return the base at the specified index. */
uint8_t Kmer::at(unsigned i) const
{
	assert(i < s_length);
	return getBaseCode(m_seq,
			seqIndexToByteNumber(i), seqIndexToBaseIndex(i));
}

/** Set the base at the specified index. */
void Kmer::set(unsigned i, uint8_t base)
{
	assert(i < s_length);
	setBaseCode(m_seq,
			seqIndexToByteNumber(i), seqIndexToBaseIndex(i), base);
}

// get a base code by the byte number and sub index
static uint8_t getBaseCode(const char* pSeq,
		unsigned byteNum, unsigned index)
{
	unsigned shiftLen = 2 * (3 - index);
	return (pSeq[byteNum] >> shiftLen) & 0x3;
}

static unsigned seqIndexToByteNumber(unsigned seqIndex)
{
	return seqIndex / 4;
}

static unsigned seqIndexToBaseIndex(unsigned seqIndex)
{
	return seqIndex % 4;
}

/** Return true if this sequence is a palindrome. */
bool Kmer::isPalindrome() const
{
	return s_length % 2 == 1 && !opt::colourSpace ? false
		: *this == ::reverseComplement(*this);
}

/** Return true if the length k-1 subsequence is a palindrome. */
bool Kmer::isPalindrome(extDirection dir) const
{
	if (s_length % 2 == 0 && !opt::colourSpace)
		return false;
	Kmer seq(*this);
	if (dir == SENSE)
		seq.shiftAppend(0);
	else
		seq.setLastBase(SENSE, 0);

	Kmer rc(*this);
	rc.reverseComplement();
	if (dir == SENSE)
		rc.setLastBase(SENSE, 0);
	else
		rc.shiftAppend(0);

	return seq == rc;
}
