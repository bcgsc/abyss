#! /usr/bin/perl

use POSIX qw(ceil floor);
use strict;

if(@ARGV != 5)
{
    die "Usage: abyss2ace.pl <kmer> <reads file> <contigs file> <align file> <outfile>\n";
}


my $LINE_LEN = 50;


my $kmer = $ARGV[0];
my $reads_file = $ARGV[1];
my $contigs_file = $ARGV[2];
my $align_file = $ARGV[3];
my $out_file = $ARGV[4];

open(OUT_ACE, ">$out_file");


#read in the contigs
my %contigs;
readFasta($contigs_file, \%contigs);

#read in the reads
my %reads;
readFasta($reads_file, \%reads);

# read in the alignments and index them by read and contig
my %readAligns;
my %contigAligns;
readAligns($align_file, \%reads, \%contigs);

#build the base segments for the contigs
#make_base_segments(\%contigs);
make_fake_base_segments(\%contigs, \%reads);


#write out the file
write_header(\%contigs, \%reads);
write_contigs(\%contigs, \%reads);
#write_reads(\%reads);

close(OUT_ACE);

#write the header line
sub write_header
{
    my($refCtgs, $refReads) = @_;
    print OUT_ACE "AS " . scalar(keys %{$refCtgs}) . " "  . scalar(keys %{$refReads}) . "\n\n";
    #print OUT_ACE "AS 1 10085\n\n";
}

sub write_contigs
{
    my($refContigs, $refReads) = @_;

    my @sKeys = sort {$a <=> $b} keys %{$refContigs};
    foreach my $ctg (@sKeys)
    {
	print "Processing contig $ctg\n";

	my $name = $ctg;
	my $l = length($refContigs->{$ctg}->{seq});

	my @aligns = @{$refContigs->{$ctg}->{aligns}};
	my @bs = @{$refContigs->{$ctg}->{bs}};


	# The alignments are guarenteed to be unique so the number of reads == number of alignments
	my $num_reads = scalar(@aligns);
	my $num_bs = scalar(@bs);

	my $oseq = chop_seq($refContigs->{$ctg}->{seq}, $LINE_LEN);

	print OUT_ACE "CO $name $l $num_reads $num_bs U\n";
	print OUT_ACE "$oseq\n\n";

	print OUT_ACE "BQ\n";

	my $qual_str = make_quality_str($l, $LINE_LEN);
	print OUT_ACE $qual_str . "\n\n";

	# write out the alignments
	# these are guarenteed to be unique to the contig
	foreach my $a (@aligns)
	{
	    my $rc_str = $a->{isRC} ? "C" : "U";
	    my $str = "AF $a->{read} $rc_str $a->{read_space_pos}";
	    print OUT_ACE "$str\n";
	}

	#write out the base segments
	my $curr_pos = 1;
	foreach my $b (@bs)
	{
	    my $end = $b->{contig_start_pos} + $b->{align_length} - 1;	
	    my $str = "BS $curr_pos $end $b->{read}";
	    print OUT_ACE "$str\n";
	    $curr_pos = $end + 1;
	}	
	print OUT_ACE "\n";

	# write out the reads
	foreach my $a (@aligns)
	{
	    my $name = $a->{read};
	    my $seq = ($a->{isRC}) ? reverseComplement($refReads->{$name}->{seq}) : $refReads->{$name}->{seq};
	    my $l = length($seq);
	    my $e = $l;
	    print OUT_ACE "RD $name $l 0 0\n";
	    print OUT_ACE chop_seq($seq, $LINE_LEN)  . "\n\n";
	    print OUT_ACE "QA 1 $e 1 $e\n";
	    print OUT_ACE "DS CHROMAT_FILE: B_RUN1_1_105_627_528 PHD_FILE: B_RUN1_1_105_627_528.phd.1 TIME: Thu Mar 19 18:38:28 1998\n\n";
	}
	#return;
    }
}

sub make_quality_str
{
    my ($num_items, $items_per_line) = @_;
    my $num_lines = ceil($num_items / $items_per_line);

    my $out = "";
    for(my $i = 0; $i < $num_lines; ++$i)
    {

	my $last_line = ($i == $num_lines - 1);

	#calculate the number of items to print
	my $num = 0;
	if($last_line)
	{
	    $num = $num_items - ($i * $items_per_line);
	}
	else
	{
	    $num = $items_per_line;
	}
	
	for(my $j = 0; $j < $num; ++$j)
	{
	    $out .= "30 ";
	}

	if(!$last_line)
	{
	    $out .= "\n";
	}
    }
    return $out;
}

sub chop_seq
{
    my ($str, $items_per_line) = @_;

    
    my $out = "";

    my $len = length($str);
    my $num_lines = ceil($len / $items_per_line);
    
    for(my $i = 0; $i < $num_lines; ++$i)
    {
	$out .= substr($str, $i*$items_per_line, $items_per_line);
	if($i != $num_lines - 1)
	{
	    $out .= "\n";
	}
    }

    return $out;
}

sub make_fake_base_segments
{
    my($refContigs, $refReads) = @_;

    print("Faking BS Seqments\n");

    my @sKeys = sort {$a <=> $b} keys %{$refContigs};
    foreach my $ctg (@sKeys)
    {
	my $full_seq = $refContigs->{$ctg}->{seq};
	for(my $i = 0; $i < length($full_seq) - $kmer + 1; ++$i)
	{
	    my $seq = substr($full_seq, $i, $kmer);
	    my $fakeRead;
	    my $fakeName = "BS_" . $ctg . "_" . $i;
	    $fakeRead->{name} = $fakeName;
	    $fakeRead->{seq} = "$seq";
	    
	    $refReads->{$fakeName} = $fakeRead;
	    
	    my $fakeAlign;
	    $fakeAlign->{read} = $fakeName;
	    $fakeAlign->{contig_start_pos} = $i + 1;
	    $fakeAlign->{read_start_pos} = 0;
	    $fakeAlign->{align_length} = $kmer;
	    $fakeAlign->{read_length} = $kmer;
	    $fakeAlign->{isRC} = 0;
	    $fakeAlign->{read_space_pos} = $fakeAlign->{contig_start_pos};
	    
	    push(@{$refContigs->{$ctg}->{aligns}}, $fakeAlign);
	    my @bs;
	    push(@bs, $fakeAlign);
	    @{$refContigs->{$ctg}->{bs}}= @bs;
	}
    }
}



#make the base segment data for the contigs by choosing a spanning set of alignments for the contig
#by definition each position in the read is covered by a k-mer so we can just choose a kmer that starts
#at each position
sub make_base_segments
{
    my($refContigs) = @_;
    
    my @sKeys = sort {$a <=> $b} keys %{$refContigs};
    foreach my $ctg (@sKeys)
    {

	print "bs contig $ctg\n";

	my @aligns = @{$refContigs->{$ctg}->{aligns}}; 

	# add the reads to each position they cover
	my $covers;
	foreach my $a (@aligns)
	{
	    my $start = $a->{contig_start_pos};
	    my $end = $start + $a->{align_length} - 1;

	    for(my $i = $start; $i <= $end; ++$i)
	    {
		#print("adding $a->{read} to " . ($i) . "\n");
		push(@{$covers->{$i}}, $a);
	    }
	}


	my $curr_pos = 1;
	my $end_pos = length($refContigs->{$ctg}->{seq}) + 1;
	my @bs;
	while($curr_pos < $end_pos)
	{
	    #print("$ctg $curr_pos $end_pos\n");
	   
	    #select a read for the current position
	    if(!defined($covers->{$curr_pos}) || scalar(@{$covers->{$curr_pos}}) <= 0)
	    {
		die "No covering read for contig $ctg at pos $curr_pos\n";
	    }
	    
	    my $r = $covers->{$curr_pos}[0];
	    #print("Selected $r->{read} $r->{contig_start_pos} $r->{align_length}\n");
	    push(@bs, $r);
	    
	    #update currpos
	    $curr_pos = $r->{contig_start_pos} + $r->{align_length};
	}

	@{$refContigs->{$ctg}->{bs}} = @bs;
	#return;
    }
}


# read the fasta file, parsing the names out
sub readFasta
{
    my ($file, $refHash) = @_;
    open(FILE, $file);
    while(my $header = <FILE>)
    {
	my $seq = <FILE>;
	
	chomp $header;
	chomp $seq;
	
	#parse the name
	my @hFields = split(' ', $header);
	my $name = substr($hFields[0], 1);
	
	#print("$name $seq\n");
	
	my $c;
	$c->{name} = $name;
	$c->{seq} = $seq;
	$refHash->{$name} = $c;
    }
    close(FILE);
}

sub readAligns
{
    my ($file, $readHash, $contigHash) = @_;
    open(FILE, $file);
    while(my $line = <FILE>)
    {
	chomp $line;

	#break the line by tabs
	my @tFields = split('\t', $line);

	my $name = $tFields[0];

	my $num_fields = @tFields;

	#consed can't handle multiple alignments
	next if($num_fields > 2);

	# Add the alignments to a temp hash by contig
	my $tempCtgHash;
	for(my $i = 1; $i < $num_fields; ++$i)
	{
	    my @aFields = split(' ', $tFields[$i]);
	    my $align;
	    $align->{read} = $name;
	    $align->{contig} = $aFields[0];

	    # ACE file positions start at 1, translate the coords
	    $align->{contig_start_pos} = $aFields[1] + 1;
	    $align->{read_start_pos} = $aFields[2];
	    $align->{align_length} = $aFields[3];
	    $align->{read_length} = $aFields[4];
	    $align->{isRC} = $aFields[5];

	    $align->{read_space_pos} = (!$align->{isRC}) ? ($align->{contig_start_pos} - $align->{read_start_pos}) :  ($align->{contig_start_pos} - (($align->{read_length} - $align->{read_start_pos}) - $align->{align_length}));

	    #index the alignments by read and contig
	    #push(@{$readHash->{$align->{read}}->{aligns}}, $align);
	    push(@{$tempCtgHash->{$align->{contig}}}, $align);
	    #push(@{$contigHash->{$align->{contig}}->{aligns}}, $align);
	}

	#Add the unique contig alignments to the main hash
	foreach my $ctg (keys %{$tempCtgHash})
	{
	    if(scalar(@{$tempCtgHash->{$ctg}}) == 1)
	    {
		push(@{$contigHash->{$ctg}->{aligns}}, $tempCtgHash->{$ctg}->[0]);
	    }
	}
    }

    close(FILE);
}

sub reverseComplement
{
    my $seq = shift;
    my @chars = split(//, $seq);
    my @comp;
    foreach my $c (@chars)
    {
        my $g = complement($c);
        push(@comp, $g);
    }

    my @rcomp = reverse(@comp);
    return join("", @rcomp);
}

sub complement
{
    my $c = shift;
    if($c eq 'A') { return 'T'; }
    if($c eq 'T') { return 'A'; }
    if($c eq 'C') { return 'G'; }
    if($c eq 'G') { return 'C'; }
    return 'N'; # shouldnt get here
}


