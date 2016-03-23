#include "BloomDBG/SpacedSeed.h"
#include <gtest/gtest.h>

using namespace std;

TEST(SpacedSeedTest, qrSeed)
{
	/*
	* Generate a Quadratic Residue (QR) seed. The background theory
	* for QR seeds is described in:
	*
	* Egidi, Lavinia, and Giovanni Manzini. "Multiple seeds
	* sensitivity using a single seed with threshold." Journal of
	* bioinformatics and computational biology 13.04 (2015): 1550011.
	*/
	ASSERT_EQ("10100011101", SpacedSeed::qrSeed(11));
}

TEST(SpacedSeedTest, qrSeedPair)
{
	/*
	 *  Generate spaced seed pattern for two mirrored QR seeds with
	 * a gap in between.
	 */
	ASSERT_EQ("101000111010000000000010111000101", SpacedSeed::qrSeedPair(33,11));
}
