#include "../Histogram.h"
#include "gtest/gtest.h"
using namespace std;

template <typename O> class MultiTypes{
	public:
	O myVar;
	void type_int() { 
		::testing::StaticAssertTypeEq<int,O>(); 
	}
	void type_size_t() { 
		::testing::StaticAssertTypeEq<size_t, O>(); 
	}
	void type_long() { 
		::testing::StaticAssertTypeEq<long, O>(); 
	}
	void type_map() {
		::testing::StaticAssertTypeEq<map<int, size_t>, O>(); 
	}
};

Histogram hi;

// test Histogram.empty()
TEST(emptyTest, base_cases) {
	// hi is just started. must be empty already.
	EXPECT_TRUE(hi.empty());
	hi.insert(2);
	EXPECT_FALSE(hi.empty());
	hi.insert(4);
	EXPECT_FALSE(hi.empty());
}

TEST(countTest, non_negative_cases){
	EXPECT_EQ(hi.size(), 2);
	hi.insert(6);
	hi.insert(8);
	hi.insert(10,5);
	EXPECT_EQ(hi.size(),9);
	EXPECT_EQ(hi.count(INT_MIN, INT_MAX),9);
	EXPECT_EQ(hi.count(8, 10),6);
	hi.insert(12);
	EXPECT_EQ(hi.size(),10);
	EXPECT_EQ(hi.count(INT_MIN, INT_MAX), 10);
	
}

TEST(sumTest, trivial_cases){
	Histogram hello;
	EXPECT_EQ(hello.sum(), 0);
}


