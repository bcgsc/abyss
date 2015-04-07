#include "Konnector/DBGBloom.h"
#include "Bloom/CascadingBloomFilter.h"
#include "Bloom/BloomFilter.h"

#include <gtest/gtest.h>
#include <string>

TEST(DBGBloom, BloomFilterPolymorphism)
{
	unsigned bits = 100000;

	Kmer::setLength(3);

	Kmer kmer1("GAC");
	 Kmer kmer2("ACC");
	  Kmer kmer3("CCA");

	CascadingBloomFilter countingBloom(bits, 2);

	countingBloom.insert(kmer1);
	countingBloom.insert(kmer1);
	countingBloom.insert(kmer2);
	countingBloom.insert(kmer2);
	countingBloom.insert(kmer3);
	countingBloom.insert(kmer3);

	DBGBloom<CascadingBloomFilter> graph(countingBloom);

	// test that expected edges exist

	boost::graph_traits< DBGBloom<CascadingBloomFilter> >::out_edge_iterator ei, ei_end;

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

	boost::graph_traits< DBGBloom<BloomFilter> >::out_edge_iterator ei2, ei_end2;

	DBGBloom<BloomFilter> graph2(countingBloom.getBloomFilter(1));

	// test that the same edges exist in non-counting
	// bloom filter

	boost::tie(ei2, ei_end2) = out_edges(kmer1, graph2);
	ASSERT_TRUE(ei2 != ei_end2);
	EXPECT_TRUE(target(*ei2, graph2) == kmer2);
	ei2++;
	EXPECT_TRUE(ei2 == ei_end2);

	boost::tie(ei2, ei_end2) = out_edges(kmer2, graph2);
	ASSERT_TRUE(ei2 != ei_end2);
	EXPECT_TRUE(target(*ei2, graph2) == kmer3);
	ei2++;
	EXPECT_TRUE(ei2 == ei_end2);
}
