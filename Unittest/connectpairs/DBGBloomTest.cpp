#include "connectpairs/DBGBloom.h"

#include <gtest/gtest.h>
#include <string>

TEST(DBGBloom, BloomFilterPolymorphism)
{
	unsigned bits = 100000;

	Kmer::setLength(3);

	Kmer kmer1("GAC");
	 Kmer kmer2("ACC");
	  Kmer kmer3("CCA");

	CountingBloomFilter countingBloom(bits);

	countingBloom.insert(kmer1);
	countingBloom.insert(kmer1);
	countingBloom.insert(kmer2);
	countingBloom.insert(kmer2);
	countingBloom.insert(kmer3);
	countingBloom.insert(kmer3);

	DBGBloom graph(countingBloom);

	// test that expected edges exist

	boost::graph_traits<DBGBloom>::out_edge_iterator ei, ei_end;

	boost::tie(ei, ei_end) = out_edges(kmer1, graph);
	ASSERT_TRUE(ei != ei_end);
	EXPECT_TRUE(target(*ei, graph) == kmer2);
	ei++;
	EXPECT_TRUE(ei == ei_end);

	boost::tie(ei, ei_end) = out_edges(kmer2, graph);
	ASSERT_TRUE(ei != ei_end);
	EXPECT_TRUE(target(*ei, graph) == kmer3);
	ei++;
	EXPECT_TRUE(ei == ei_end);

	DBGBloom graph2(countingBloom.getBloomFilter(1));

	// test that the same edges exist in non-counting
	// bloom filter

	boost::tie(ei, ei_end) = out_edges(kmer1, graph2);
	ASSERT_TRUE(ei != ei_end);
	EXPECT_TRUE(target(*ei, graph2) == kmer2);
	ei++;
	EXPECT_TRUE(ei == ei_end);

	boost::tie(ei, ei_end) = out_edges(kmer2, graph2);
	ASSERT_TRUE(ei != ei_end);
	EXPECT_TRUE(target(*ei, graph2) == kmer3);
	ei++;
	EXPECT_TRUE(ei == ei_end);
}
