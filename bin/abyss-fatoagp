#!/usr/bin/env perl
# Convert a FASTA file of scaffolds to a FASTA file of contigs and an
# AGP file.
# Written by Shaun Jackman <sjackman@bcgsc.ca>.

use strict;
use Getopt::Std qw'getopts';

my %opt;
getopts 'f:s:S:', \%opt;
my $opt_fasta = $opt{'f'};
# scaffolds shorter than this length will be excluded
my $opt_min_scaf_len = defined $opt{'s'} ? $opt{'s'} : 200;
# scaftigs shorter than this length will be masked with "N"s
my $opt_min_ctg_len = defined $opt{'S'} ? $opt{'S'} : 50;

open FASTA, ">$opt_fasta"
	or die "error: `$opt_fasta': $!\n"
	if $opt_fasta;

while (<>) {
	die unless /^>/;
	chomp;
	my ($scafid, undef) = split ' ', $_, 2;
	substr $scafid, 0, 1, '';

	my $scafseq = <>;
	chomp $scafseq;
	my $scaflen = $scafseq =~ tr/ACGTacgt//;
	next if $scaflen < $opt_min_scaf_len;

	# mask scaftigs shorter than length threshold with "N"s
	my @ctgseqs = split /([Nn]+)/, $scafseq;
	foreach my $ctgseq (@ctgseqs) {
		next if /^[nN]/;
		if (length($ctgseq) < $opt_min_ctg_len) {
			$ctgseq = "N" x length($ctgseq);
		}
	}
	# rejoin and split to merge adjacent stretches of "N"s
	$scafseq = join '', @ctgseqs;
	next unless $scafseq =~ /[^nN]/;
	@ctgseqs = split /([Nn]+)/, $scafseq;

	my $i = 0;
	my $x = 0;
	for my $ctgseq (@ctgseqs) {
		my $len = length $ctgseq;
		$i++ if ($len == 0);
		next if ($len == 0);
		# object object_beg object_end part_number
		print 'scaffold', $scafid, "\t",
			$x + 1, "\t",
			$x + $len, "\t",
			$i + 1, "\t";
		if ($ctgseq =~ /^[nN]/) {
			# component_type gap_length gap_type linkage
			print "N\t", $len, "\tscaffold\tyes\tpaired-ends\n";
		} else {
			my $ctgid = 'contig' . $scafid . '_' . ($i / 2);
			# component_type component_id
			# component_beg component_end orientation
			print "W\t", $ctgid, "\t1\t", $len, "\t+\n";
			print FASTA '>', $ctgid, "\n", $ctgseq, "\n"
				if $opt_fasta;
		}
		$i++;
		$x += $len;
	}
}
