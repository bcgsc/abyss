#include "Common/KmerIterator.h"
#include <gtest/gtest.h>

TEST(KmerIteratorTest, NoIllegalChars)
{
	unsigned k = 3;
	Kmer::setLength(k);
	KmerIterator i("AGCTA", k);

	ASSERT_EQ(Kmer("AGC"), *i);
	i++;
	ASSERT_EQ(Kmer("GCT"), *i);
	i++;
	ASSERT_EQ(Kmer("CTA"), *i);
	i++;
	ASSERT_EQ(KmerIterator::end(), i);
}

TEST(KmerIteratorTest, IllegalChars)
{
	unsigned k = 3;
	Kmer::setLength(k);
	KmerIterator i("AGCTNTAG", k);

	ASSERT_EQ(Kmer("AGC"), *i);
	i++;
	ASSERT_EQ(Kmer("GCT"), *i);
	i++;
	ASSERT_EQ(Kmer("TAG"), *i);
	i++;
	ASSERT_EQ(KmerIterator::end(), i);
}

TEST(KmerIteratorTest, SeqLengthLessThanK)
{
	unsigned k = 3;
	Kmer::setLength(k);
	KmerIterator i("AG", k);
	ASSERT_EQ(KmerIterator::end(), i);
}
