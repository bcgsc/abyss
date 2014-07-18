#include "Common/BitUtil.h"
#include "gtest/gtest.h"
#include <boost/utility/binary.hpp>

/** Test limits */
TEST(popcountTest, boundaries)
{
	EXPECT_EQ(64ULL, popcount(0xffffffffffffffffULL));
	EXPECT_EQ(0ULL, popcount(0ULL));
}

/** Test some random values */
TEST(popcountTest, random_values)
{
	EXPECT_EQ(45ULL, popcount(0x992E54FFFFFFFBA1ULL));
	EXPECT_EQ(45ULL, popcount(0x814BC5FFFFFFF7FULL));
	EXPECT_EQ(46ULL, popcount(0x815BC5FFFFFFF7FULL));
}

TEST(readBitsTest, overwriteBits)
{
	// NOTE!: In the binary diagrams below, the rightmost bit of each byte
	// is the LSB.

	size_t bitSize, bitOffset;
	char src[2];
	char dest[4];

	// SUBTEST: src size < 1 byte
	//
	// src             = 000
	// bit offset      = 9
	// dest            = 11111111 11111111 11111111 11111111
	// expected result = 11111111 10001111 11111111 11111111

	bitSize = 3;
	bitOffset = 9;
	memset(src, 0x00, 2);
	memset(dest, 0xFF, 4);

	std::istringstream in1(std::string(src, (bitSize + 7)/8));

	readBits(in1, dest, bitSize, bitOffset);

	EXPECT_EQ((char)0xFF, dest[0]);
	EXPECT_EQ((char)BOOST_BINARY(10001111), dest[1]);
	EXPECT_EQ((char)0xFF, dest[2]);
	EXPECT_EQ((char)0xFF, dest[3]);

	// SUBTEST: byte-aligned with partial last byte
	//
	// src             = 00000000 000
	// bit offset      = 8
	// dest            = 11111111 11111111 11111111 11111111
	// expected result = 11111111 00000000 00011111 11111111

	bitSize = 11;
	bitOffset = 8;
	memset(src, 0x00, 2);
	memset(dest, 0xFF, 4);

	std::istringstream in2(std::string(src, (bitSize + 7)/8));

	readBits(in2, dest, bitSize, bitOffset);

	EXPECT_EQ((char)0xFF, dest[0]);
	EXPECT_EQ((char)0x00, dest[1]);
	EXPECT_EQ((char)BOOST_BINARY(00011111), dest[2]);
	EXPECT_EQ((char)0xFF, dest[3]);

	// SUBTEST: not byte-aligned, src size > 1 byte
	//
	// src             = 00000000 00000000
	// bit offset      = 4
	// dest            = 11111111 11111111 11111111 11111111
	// expected result = 11110000 00000000 00011111 11111111

	bitSize = 15;
	bitOffset = 4;
	memset(src, 0x00, 2);
	memset(dest, 0xFF, 4);

	std::istringstream in3(std::string(src, (bitSize + 7)/8));

	readBits(in3, dest, bitSize, bitOffset);

	EXPECT_EQ((char)BOOST_BINARY(11110000), dest[0]);
	EXPECT_EQ((char)0x00, dest[1]);
	EXPECT_EQ((char)BOOST_BINARY(00011111), dest[2]);
	EXPECT_EQ((char)0xFF, dest[3]);
}

TEST(readBitsTest, orBits)
{
	// NOTE!: In the binary diagrams below, the rightmost bit of each byte
	// is the LSB.

	size_t bitSize, bitOffset;
	char src[2];
	char dest[4];

	// SUBTEST: src size < 1 byte
	//
	// src             = 101
	// bit offset      = 9
	// dest            = 10101010 10101010 10101010 10101010
	// expected result = 10101010 11111010 10101010 10101010

	bitSize = 3;
	bitOffset = 9;
	memset(src, 0x00, 2);
	src[0] = BOOST_BINARY(10100000);
	memset(dest, BOOST_BINARY(10101010), 4);

	std::istringstream in1(std::string(src, (bitSize + 7)/8));

	readBits(in1, dest, bitSize, bitOffset, BITWISE_OR);

	EXPECT_EQ((char)BOOST_BINARY(10101010), dest[0]);
	EXPECT_EQ((char)BOOST_BINARY(11111010), dest[1]);
	EXPECT_EQ((char)BOOST_BINARY(10101010), dest[2]);
	EXPECT_EQ((char)BOOST_BINARY(10101010), dest[3]);

	// SUBTEST: byte-aligned with partial last byte
	//
	// src             = 01010101 010
	// bit offset      = 8
	// dest            = 10101010 10101010 10101010 10101010
	// expected result = 10101010 11111111 11101010 10101010

	bitSize = 11;
	bitOffset = 8;
	memset(src, 0x00, 2);
	src[0] = BOOST_BINARY(01010101);
	src[1] = BOOST_BINARY(01000000);
	memset(dest, BOOST_BINARY(10101010), 4);

	std::istringstream in2(std::string(src, (bitSize + 7)/8));

	readBits(in2, dest, bitSize, bitOffset, BITWISE_OR);

	EXPECT_EQ((char)BOOST_BINARY(10101010), dest[0]);
	EXPECT_EQ((char)BOOST_BINARY(11111111), dest[1]);
	EXPECT_EQ((char)BOOST_BINARY(11101010), dest[2]);
	EXPECT_EQ((char)BOOST_BINARY(10101010), dest[3]);

	// SUBTEST: not byte-aligned, src size > 1 byte
	//
	// src             = 01010101 0101010
	// bit offset      = 4
	// dest            = 10101010 10101010 10101010 10101010
	// expected result = 10101111 11111111 11101010 10101010

	bitSize = 15;
	bitOffset = 4;
	memset(src, 0x00, 2);
	src[0] = BOOST_BINARY(01010101);
	src[1] = BOOST_BINARY(01010100);
	memset(dest, BOOST_BINARY(10101010), 4);

	std::istringstream in3(std::string(src, (bitSize + 7)/8));

	readBits(in3, dest, bitSize, bitOffset, BITWISE_OR);

	EXPECT_EQ((char)BOOST_BINARY(10101111), dest[0]);
	EXPECT_EQ((char)BOOST_BINARY(11111111), dest[1]);
	EXPECT_EQ((char)BOOST_BINARY(11101010), dest[2]);
	EXPECT_EQ((char)BOOST_BINARY(10101010), dest[3]);
}

TEST(readBitsTest, andBits)
{
	// NOTE!: In the binary diagrams below, the rightmost bit of each byte
	// is the LSB.

	size_t bitSize, bitOffset;
	char src[2];
	char dest[4];

	// SUBTEST: src size < 1 byte
	//
	// src             = 101
	// bit offset      = 9
	// dest            = 11111111 11111111 11111111 11111111
	// expected result = 11111111 11011111 11111111 11111111

	bitSize = 3;
	bitOffset = 9;
	memset(src, 0x00, 2);
	src[0] = BOOST_BINARY(10100000);
	memset(dest, BOOST_BINARY(11111111), 4);

	std::istringstream in1(std::string(src, (bitSize + 7)/8));

	readBits(in1, dest, bitSize, bitOffset, BITWISE_AND);

	EXPECT_EQ((char)BOOST_BINARY(11111111), dest[0]);
	EXPECT_EQ((char)BOOST_BINARY(11011111), dest[1]);
	EXPECT_EQ((char)BOOST_BINARY(11111111), dest[2]);
	EXPECT_EQ((char)BOOST_BINARY(11111111), dest[3]);

	// SUBTEST: byte-aligned with partial last byte
	//
	// src             = 01010101 010
	// bit offset      = 8
	// dest            = 11111111 11111111 11111111 11111111
	// expected result = 11111111 01010101 01011111 11111111

	bitSize = 11;
	bitOffset = 8;
	memset(src, 0x00, 2);
	src[0] = BOOST_BINARY(01010101);
	src[1] = BOOST_BINARY(01000000);
	memset(dest, BOOST_BINARY(11111111), 4);

	std::istringstream in2(std::string(src, (bitSize + 7)/8));

	readBits(in2, dest, bitSize, bitOffset, BITWISE_AND);

	EXPECT_EQ((char)BOOST_BINARY(11111111), dest[0]);
	EXPECT_EQ((char)BOOST_BINARY(01010101), dest[1]);
	EXPECT_EQ((char)BOOST_BINARY(01011111), dest[2]);
	EXPECT_EQ((char)BOOST_BINARY(11111111), dest[3]);

	// SUBTEST: not byte-aligned, src size > 1 byte
	//
	// src             = 01010101 0101010
	// bit offset      = 4
	// dest            = 11111111 11111111 11111111 11111111
	// expected result = 11110101 01010101 01011111 11111111

	bitSize = 15;
	bitOffset = 4;
	memset(src, 0x00, 2);
	src[0] = BOOST_BINARY(01010101);
	src[1] = BOOST_BINARY(01010100);
	memset(dest, BOOST_BINARY(11111111), 4);

	std::istringstream in3(std::string(src, (bitSize + 7)/8));

	readBits(in3, dest, bitSize, bitOffset, BITWISE_AND);

	EXPECT_EQ((char)BOOST_BINARY(11110101), dest[0]);
	EXPECT_EQ((char)BOOST_BINARY(01010101), dest[1]);
	EXPECT_EQ((char)BOOST_BINARY(01011111), dest[2]);
	EXPECT_EQ((char)BOOST_BINARY(11111111), dest[3]);
}
