#include "Common/StringUtil.h"
#include "gtest/gtest.h"
using namespace std;

TEST(chop_test, base_case_2)
{
	string myString = "me";
	EXPECT_EQ((unsigned)2, myString.length());
	EXPECT_EQ('e', chop(myString));
	EXPECT_EQ((unsigned)1, myString.length());
}

TEST(chop_test, trivial_gt_length2)
{
	string myString = "something";
	EXPECT_EQ((unsigned)9, myString.length());
	EXPECT_EQ('g', chop(myString));
	EXPECT_EQ((unsigned)8, myString.length());
	EXPECT_EQ('n', chop(myString));
	EXPECT_EQ((unsigned)7, myString.length());
}

TEST(chomp_test, base_cases)
{
	// test .length=1
	string myString = "a";
	EXPECT_TRUE(1 == myString.length());
	EXPECT_FALSE(chomp(myString));
	string anotherString = "\n";
	EXPECT_TRUE(1 == anotherString.length());
	EXPECT_TRUE(chomp(anotherString));

	// test .length=2
	string greatString = "ab";
	EXPECT_FALSE(chomp(greatString));
	EXPECT_EQ((unsigned)2, greatString.length());

	string badString = "a\n";
	EXPECT_TRUE(chomp(badString));
	EXPECT_EQ((unsigned)1, badString.length());
}

TEST(toSI_test, all_the_cases)
{
	// negative and zero values
	EXPECT_EQ ("-0.000123 ", toSI(-0.0001234));
	EXPECT_EQ("1e-13 ", toSI(0.0000000000001));
	EXPECT_EQ("-1.2e-13 ", toSI(-0.00000000000012));
	EXPECT_EQ("-1.23e-13 ", toSI(-0.000000000000123456));
	EXPECT_EQ("0 ", toSI(0));
	EXPECT_EQ("-0 ", toSI(-0.000));

	//trivially all posible values
	EXPECT_EQ("1.23 k", toSI(1234));
	EXPECT_EQ("123 M", toSI(123440111));
	EXPECT_EQ("23.4 M", toSI(23440111));
	EXPECT_EQ("1.23 G", toSI(1234123123));
	EXPECT_EQ("123 G", toSI(123440111222));
	EXPECT_EQ("23.4 G", toSI(23440222111));
}

TEST(fromSI_test, all_the_cases)
{
	// zero values
	EXPECT_EQ(0, fromSI("0"));
	EXPECT_EQ(0, fromSI("0T"));

	// negative values
	EXPECT_EQ(-1.23e3, fromSI("-1.23k"));
	EXPECT_EQ(-9.06e6, fromSI("-9.06M"));
	EXPECT_EQ(-1.234e9, fromSI("-1.234G"));
	EXPECT_EQ(-1.234567e12, fromSI("-1.234567T"));

	// positive values
	EXPECT_EQ(1.23e3, fromSI("1.23k"));
	EXPECT_EQ(9.06e6, fromSI("9.06M"));
	EXPECT_EQ(1.234e9, fromSI("1.234G"));
	EXPECT_EQ(1.234567e12, fromSI("1.234567T"));
}

TEST(bytesToSI_test, all_the_cases)
{
	// zero values
	EXPECT_EQ("0", bytesToSI(0));

	// unit conversion
	EXPECT_EQ("1", bytesToSI(1));
	EXPECT_EQ("1k", bytesToSI(1024));
	EXPECT_EQ("1M", bytesToSI(1048576UL));
	EXPECT_EQ("1G", bytesToSI(1073741824UL));
	EXPECT_EQ("1T", bytesToSI(1099511627776ULL));
}

template <typename T>
class MultiTypes{
 public:
 	T myVar;
	void type_int() { ::testing::StaticAssertTypeEq<int, T>(); }
	void type_double() { ::testing::StaticAssertTypeEq<double, T>(); }
	void type_string() { ::testing::StaticAssertTypeEq<string, T>(); }
};

TEST(toEng_test, integer_cases)
{
	EXPECT_EQ("1234", toEng(1234));

	MultiTypes<int> temp;
	temp.type_int();
	temp.myVar = 1234;
	EXPECT_EQ(temp.myVar, 1234);
	EXPECT_EQ("1234", toEng(temp.myVar));
	temp.myVar = 12345678;

	EXPECT_EQ("12.35e6", toEng(temp.myVar));
	temp.myVar = 123456789;
	EXPECT_EQ ("123.5e6" , toEng(temp.myVar));
}

TEST(toEng_test, double_cases)
{
	MultiTypes<double> temp;
	temp.type_double();

	temp.myVar = 123.456;
	EXPECT_EQ (temp.myVar, 123.456);
	EXPECT_EQ ("123.5", toEng(temp.myVar));
	temp.myVar = 123456789.9;
	EXPECT_EQ ("123.5e6", toEng(temp.myVar));
}

TEST(toEng_test, string_cases)
{
	MultiTypes<string> temp;
	temp.type_string();
	temp.myVar = "123.456";
}

TEST (startsWith_test, trivial_cases)
{
	EXPECT_TRUE (startsWith("hello", "hell"));
	EXPECT_TRUE (startsWith("hello", "h"));
	EXPECT_FALSE (startsWith("hello", "hello"));
	EXPECT_TRUE (startsWith("hello", ""));
	EXPECT_FALSE (startsWith ("whatever", "who"));
}

TEST(endsWith_test, any_cases)
{
	// suffix should not be the string itself
	EXPECT_FALSE(endsWith("hello", "hello"));
	EXPECT_FALSE(endsWith("", ""));
	EXPECT_TRUE(endsWith("hello", ""));

	// EXPECT_FALSE(endsWith("hello", NULL));
	// NULL is not valid
	EXPECT_TRUE(endsWith("hello", "ello"));
	EXPECT_TRUE(endsWith("hello", "o"));
	EXPECT_FALSE(endsWith("hell", "hello"));
	EXPECT_FALSE(endsWith("hell", "heaven"));
}

TEST(SIToBytes_test, unit_conversions)
{
	EXPECT_EQ(1024u, SIToBytes("1024"));
	EXPECT_EQ(1024u, SIToBytes("1k"));
	EXPECT_EQ(1536u, SIToBytes("1.5k"));
	EXPECT_EQ(1048576u, SIToBytes("1M"));
	EXPECT_EQ(1073741824u, SIToBytes("1G"));
	EXPECT_EQ(1099511627776u, SIToBytes("1T"));
}

TEST(SIToBytes_test, error_handling)
{
	unsigned long long bytes;

	// non-number

	istringstream nonNumber("not-a-number");
	bytes = SIToBytes(nonNumber);
	EXPECT_EQ(0u, bytes);
	EXPECT_TRUE(nonNumber.fail());

	// unrecognized suffix

	istringstream illegalUnits("1024y");
	bytes = SIToBytes(nonNumber);
	EXPECT_EQ(0u, bytes);
	EXPECT_TRUE(nonNumber.fail());

	// valid string should set eof

	istringstream valid("500M");
	bytes = SIToBytes(valid);
	EXPECT_EQ(524288000u, bytes);
	EXPECT_TRUE(valid.eof());
}
