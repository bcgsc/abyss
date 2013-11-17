#include "Common/BitUtil.h"
#include "gtest/gtest.h"

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
