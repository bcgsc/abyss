#include <sstream>
#include <string>

#include "Common/SAM.h"
#include "gtest/gtest.h"

/** Verify SAM records are stored and used correctly */

using namespace std;

/** @Return whether two Alignments are equivalent. */
bool operator==(const Alignment& a, const Alignment& b)
{
	return a.read_start_pos == b.read_start_pos &&
		a.align_length == b.align_length &&
		a.read_length == b.read_length;
}


// Test SAM::parseCigar()
TEST(parseCigar, check_alignment)
{
	Alignment a;
	a.align_length = 40;
	a.read_start_pos = 0;
	a.read_length = 40;
	EXPECT_EQ(a, SAMAlignment::parseCigar("40M", false));

	a.align_length = 40;
	a.read_start_pos = 20;
	a.read_length = 60;
	EXPECT_EQ(a, SAMAlignment::parseCigar("20S40M", false));

	a.align_length = 40;
	a.read_start_pos = 0;
	a.read_length = 60;
	EXPECT_EQ(a, SAMAlignment::parseCigar("40M20S", false));
	EXPECT_EQ(a, SAMAlignment::parseCigar("20S40M", true));

	a.align_length = 40;
	a.read_start_pos = 20;
	a.read_length = 70;
	EXPECT_EQ(a, SAMAlignment::parseCigar("20S40M10S", false));
	EXPECT_EQ(a, SAMAlignment::parseCigar("10S40M20S", true));

	a.align_length = 40;
	a.read_start_pos = 20;
	a.read_length = 70;
	EXPECT_EQ(a, SAMAlignment::parseCigar("20I40M10S", false));

	a.align_length = 40;
	a.read_start_pos = 30;
	a.read_length = 80;
	EXPECT_EQ(a, SAMAlignment::parseCigar("20M10I40M10S", false));

	a.align_length = 40;
	a.read_start_pos = 0;
	a.read_length = 80;
	EXPECT_EQ(a, SAMAlignment::parseCigar("40M10I20M10S", false));
}

// Check that we error when an invalid CIGAR is given.
TEST(parseCigarDeath, invalid_cigar)
{
	EXPECT_DEATH(SAMAlignment::parseCigar("20SS", false), "error: invalid CIGAR: `20SS'");
	EXPECT_DEATH(SAMAlignment::parseCigar("20m", false), "error: invalid CIGAR: `20m'");
}

// Test friend std::istream& operator >>(std::istream& in, SAMRecord& o)
TEST(parseSAMInput, check_values)
{
	std::string testSAMString(  "1:497:R:-272+13M17D24M	113	1	497	37	37M	15	100338662	0	"
								"CGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAG	0;==-==9;>>>>>=>>>>>>>>>>>=>>>>>>>>>>	"
								"XT:A:U	NM:i:0	SM:i:37	AM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:37"
								"\n" 
								"19:20389:F:275+18M2D19M	99	1	17644	0	37M	=	17919	314	"
								"TATGACTGCTAATAATACCTACACATGTTAGAACCAT	>>>>>>>>>>>>>>>>>>>><<>>><<>>4::>>:<9	"
								"RG:Z:UM0098:1	XT:A:R	NM:i:0	SM:i:0	AM:i:0	X0:i:4	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:37"
								);
	std::istringstream testSAMIStringStream(testSAMString);
	SAMRecord testSAMRecord;
	testSAMIStringStream >> testSAMRecord;
	EXPECT_EQ("1:497:R:-272+13M17D24M", testSAMRecord.qname);
	EXPECT_EQ(113, testSAMRecord.flag);
	EXPECT_EQ("1", testSAMRecord.rname);
	EXPECT_EQ(496, testSAMRecord.pos);
	EXPECT_EQ(37, testSAMRecord.mapq);
	EXPECT_EQ("37M", testSAMRecord.cigar);
	EXPECT_EQ("15", testSAMRecord.mrnm);
	EXPECT_EQ(100338661, testSAMRecord.mpos);
	EXPECT_EQ(0, testSAMRecord.isize);
#if SAM_SEQ_QUAL
	EXPECT_EQ("CGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAG", testSAMRecord.seq);
	EXPECT_EQ("0;==-==9;>>>>>=>>>>>>>>>>>=>>>>>>>>>>", testSAMRecord.qual);
	EXPECT_EQ("XT:A:U	NM:i:0	SM:i:37	AM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:37", testSAMRecord.tags);
#endif
	testSAMIStringStream >> testSAMRecord;
	EXPECT_EQ("19:20389:F:275+18M2D19M", testSAMRecord.qname);
	EXPECT_EQ(99, testSAMRecord.flag);
	EXPECT_EQ("1", testSAMRecord.rname);
	EXPECT_EQ(17643, testSAMRecord.pos);
	EXPECT_EQ(0, testSAMRecord.mapq);
	EXPECT_EQ("37M", testSAMRecord.cigar);
	EXPECT_EQ("1", testSAMRecord.mrnm);
	EXPECT_EQ(17918, testSAMRecord.mpos);
	EXPECT_EQ(314, testSAMRecord.isize);
#if SAM_SEQ_QUAL
	EXPECT_EQ("TATGACTGCTAATAATACCTACACATGTTAGAACCAT", testSAMRecord.seq);
	EXPECT_EQ(">>>>>>>>>>>>>>>>>>>><<>>><<>>4::>>:<9", testSAMRecord.qual);
	EXPECT_EQ("RG:Z:UM0098:1	XT:A:R	NM:i:0	SM:i:0	AM:i:0	X0:i:4	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:37", testSAMRecord.tags);
#endif
	
}
