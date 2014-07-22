#include "Common/SAM.h"
#include <gtest/gtest.h>

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
