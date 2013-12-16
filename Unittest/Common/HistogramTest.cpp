#include "Common/Histogram.h"
#include "gtest/gtest.h"


// test Histogram.empty()
TEST(emptyTest, base_cases)
{
	Histogram hi;
	EXPECT_TRUE(hi.empty());
	hi.insert(2);
	EXPECT_FALSE(hi.empty());
	hi.insert(4);
	EXPECT_FALSE(hi.empty());
}

// test Histogram.count()
TEST(countTest, non_negative_cases)
{
	Histogram hi;
	hi.insert(2);
	hi.insert(4);
	EXPECT_EQ(hi.size(), (unsigned)2);
	hi.insert(6);
	hi.insert(8);
	hi.insert(10, 5);
	EXPECT_EQ(hi.size(), (unsigned)9);
	EXPECT_EQ(hi.count(INT_MIN, INT_MAX), (unsigned)9);
	EXPECT_EQ(hi.count(8, 10), (unsigned)6);
	hi.insert(12);
	EXPECT_EQ(hi.size(), (unsigned)10);
	EXPECT_EQ(hi.count(INT_MIN, INT_MAX), (unsigned)10);
}

// test Histogram.sum()
TEST(sumTest, trivial_cases)
{
	Histogram hello;
	EXPECT_EQ(hello.sum(), (unsigned)0);
}
