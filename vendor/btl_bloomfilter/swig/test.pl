#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin/./";
use BloomFilter;

#Same inputs as the adhoc tests
$filterSize = 1000000000;

#Check if filter is able to report expected results
$filter = BloomFilter::BloomFilter->new($filterSize, 5, 20);

$filter->insert("ATCGGGTCATCAACCAATAT");
$filter->insert("ATCGGGTCATCAACCAATAC");
$filter->insert("ATCGGGTCATCAACCAATAG");
$filter->insert("ATCGGGTCATCAACCAATAA");

if (!$filter->contains("ATCGGGTCATCAACCAATAT")
	&&!$filter->contains("ATCGGGTCATCAACCAATAC")
	&& !$filter->contains("ATCGGGTCATCAACCAATAG")
	&& !BloomFilter::BloomFilter::contains($filter, "ATCGGGTCATCAACCAATAA")) {
	print "Filter did not contain expected. \n";
}

if ($filter->contains("ATCGGGTCATCAACCAATTA")
	&& $filter->("ATCGGGTCATCAACCAATTC")) {
	print "Filter contained unexpected. \n";
}

print "de novo bf tests done \n";

#Check storage can occur properly
$fileName = "BloomFilter.bf";
$filter->storeFilter($fileName);

$filter2 = new BloomFilter::BloomFilter($fileName);

if (!$filter2->contains("ATCGGGTCATCAACCAATAT")
	&&!$filter2->contains("ATCGGGTCATCAACCAATAC")
	&& !$filter2->contains("ATCGGGTCATCAACCAATAG")
	&& !$filter2->contains("ATCGGGTCATCAACCAATAA")) {
	print "Filter2 did not contain expected. \n";
}

if ($filter2->contains("ATCGGGTCATCAACCAATTA")
	&& $filter2->contains("ATCGGGTCATCAACCAATTC")) {
	print "Filter2 contained unexpected. \n";

}
print "premade bf tests done\n";

$pop = $filter2->getPop();
$hash = $filter2->getHashNum();
$ksize = $filter2->getKmerSize();
$size = $filter2->getFilterSize();
print "Filter Info: Pop - $pop, numHash - $hash, kmerSize - $ksize, size - $size\n";

#RollingHashIterator tests
my $k = 5;
$str = "TAGAATCACCCAAAGA";
$bloom = new BloomFilter::BloomFilter(10000, 4, $k);
#$itr = new BloomFilter::RollingHashIterator($str, 4, $k);

BloomFilter::insertSeq($bloom, $str, 4, $k);

#my $count = 0;
#my $next = $itr->getNext();
#while ($next) {
#	$bloom->insert($next);
#    print substr($str, $count, $k) . " " .  $count++ . "\n";
#	$next =  $itr->getNext();
#}

for (my $i = 0; $i < length($str) - $k + 1; $i++) {
	my $kmer = substr($str, $i, $k);
	print $i . " ";
	if($bloom->contains($kmer)){
		print $kmer . " found\n";
	}
	else{
		print $kmer . " not found\n";
	}
}

print "Done!\n";

exit;
