#!/usr/bin/env perl
# abyss-fac (FASTA count)
# Calculate assembly contiguity statistics, such as N50.
# Written by Shaun Jackman <sjackman@bcgsc.ca>.
use strict;
use Getopt::Std qw'getopts';

$| = 1;

my %opt;
getopts 'g:hHjt:', \%opt;
my $opt_threshold = defined $opt{'t'} ? $opt{'t'} : 200;
my $opt_filename = $opt{'H'} || (@ARGV > 1 && !$opt{'h'});
my $opt_jira = $opt{'j'};
my $opt_genome_size = $opt{'g'};

sub eng($)
{
	my $x = shift;
	return $x if $x < 10000000;
	return substr($x / 1000000, 0, 5) . 'e6' if $x < 1000000000;
	return substr($x / 1000000000, 0, 5) . 'e9';
}

my ($short, $sum);
my @x;

sub count($$)
{
	my $id = shift;
	my $seq = uc shift;
	my $x = $seq =~ tr/ACGT//;
	my $colourspace = $seq =~ tr/0123//;
	die unless $x == 0 || $colourspace == 0;
	$x = $colourspace if $x == 0;
	if ($x < $opt_threshold) {
		$short++;
		return;
	}
	$sum += $x;
	push @x, $x;
}

sub fac($)
{
	my $path = shift;
	$short = $sum = 0;
	@x = ();

	my $id;
	my $seq;
	open IN, "<$path" or die "$path: $!\n";
	while (<IN>) {
		chomp;
		if (/^>/) {
			count $id, $seq if defined $id;
			$id = $_;
			$seq = '';
		} else {
			$seq .= $_;
		}
	}
	count $id, $seq if defined $id;
	close IN;

	my $n = @x;
	if ($n > 0) {
		@x = sort { $a <=> $b } @x;
		my $min = $x[0];
		my $max = $x[-1];

		my $n50_target = defined $opt_genome_size
				? $opt_genome_size : $sum;
		my ($n20, $n20sum,
			$nn50, $n50, $n50sum,
			$n80, $n80sum);
		while (@x > 0 && $n80sum < 0.8 * $n50_target) {
			my $x = pop @x;
			if ($n20sum < 0.2 * $n50_target) {
				$n20 = $x;
				$n20sum += $x;
			}
			if ($n50sum < 0.5 * $n50_target) {
				$nn50++;
				$n50 = $x;
				$n50sum += $x;
			}
			if ($n80sum < 0.8 * $n50_target) {
				$n80 = $x;
				$n80sum += $x;
			}
		}

		my $ntotal = $short + $n;
		format Spaces =
@<<<<<<<@<<<<<<<@<<<<<<<@<<<<<<<@<<<<<<<@<<<<<<<@<<<<<<<@<<<<<<<@<<<<<<<@*
eng($ntotal), eng($n), $nn50, $min, $n80, $n50, $n20, $max, eng($sum), $path
.
		format Pipes =
|@<<<<<<|@<<<<<<|@<<<<<<|@<<<<<<|@<<<<<<|@<<<<<<|@<<<<<<|@<<<<<<|@<<<<<<|@*|
eng($ntotal), eng($n), $nn50, $min, $n80, $n50, $n20, $max, eng($sum), $path
.
		$~ = $opt_jira ? 'Pipes' : 'Spaces';
		$^ = $opt_jira ? 'Pipes_TOP' : 'Spaces_TOP';
		write;
	} else {
		print STDERR "warning: `$path' is empty\n";
	}
}

format Spaces_TOP =
n       n:@<<<< L50   min     N80     N50     N20     max     sum
$opt_threshold
.
format Pipes_TOP =
||n    ||n:@<<<||L50 ||min   ||N80   ||N50   ||N20   ||max   ||sum   ||
$opt_threshold
.

@ARGV = ('-') if @ARGV == 0;
fac $_ foreach @ARGV;
