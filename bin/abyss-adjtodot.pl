#!/usr/bin/env perl
# Convert an ABySS adjacency file to GraphViz dot format.
# Written by Shaun Jackman <sjackman@bcgsc.ca>.
use strict;

print "digraph adj {\n";

while (<>) {
	chomp;
	my ($id, $length, $coverage, $a, $b);
	if (/;.*;/) {
		($id, $length, $coverage, $a, $b)
			= /^([^ ]+)\s+([^ ]+)\s+([^ ]+)\s;\s*(.*)\s;\s*(.*)$/;
	} elsif (/;/) {
		($id, $length, $a, $b)
			= /^([^ ]+)\s+([^ ]+)\s*(.*)\s;\s*(.*)$/;
	} else {
		s/,0/+/g;
		s/,1/-/g;
		($id, $length, $a, $b) = /(.*) (.*) \[(.*)\] \[(.*)\]/;
	}
	my @a = split ' ', $a;
	my @b = split ' ', $b;

	my $attr = "l=$length";
	$attr .= " C=$coverage" if defined $coverage;
	print qq{"$id+" \[$attr];\n};

	print qq{"$id+"};
	if (@a > 0) {
		print ' -> {';
		print qq{ "$_"} for @a;
		print ' }';
	}
	print ";\n";

	print qq{"$id-" \[$attr];\n};

	print qq{"$id-"};
	if (@b > 0) {
		print ' -> {';
		for (@b) {
			my $x = $_;
			$x =~ y/+-/-+/;
			print qq{ "$x"};
		}
		print ' }';
	}
	print ";\n";
}

print "}\n";
