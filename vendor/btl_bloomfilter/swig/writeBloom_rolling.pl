#!/usr/bin/env perl


#AUTHOR
#   Rene Warren
#   rwarren at bcgsc.ca

#NAME
#writeBloom

#SYNOPSIS

#DOCUMENTATION
#   LINKS-readme.txt distributed with this software @ www.bcgsc.ca
#   http://www.bcgsc.ca/platform/bioinfo/software/links
#   We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca
#   If you use LINKS, the LINKS code or ideas, please cite our work

#LICENSE
#   LINKS Copyright (c) 2014-2015 Canada's Michael Smith Genome Science Centre.  All rights reserved.

#LINKS is released under the BC Cancer Agency software license agreement (academic use). Details of the license can be accessed at: http://www.bcgsc.ca/platform/bioinfo/license/bcca_2010
#For commercial use, please contact rwarren@bcgsc.ca


use strict;
use POSIX;
use FindBin;
#perl 5.8.8
#use lib "$FindBin::Bin/./lib/Bloom-Faster-1.6/bloom/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi";
#perl 5.10.0
#use lib "$FindBin::Bin/../lib/Bloom-Faster-1.7/bloom5-16-3/lib/site_perl/5.16.3/x86_64-linux";###rebuild against PERL version or replace by pre-built
#use Bloom::Faster;
use lib "$FindBin::Bin/./";
use BloomFilter;
use Getopt::Std;
use Net::SMTP;
use Time::HiRes;
use vars qw($opt_f $opt_k $opt_p);
getopts('f:k:p:');
my ($k,$fpr)=(15,0.0001);

#-------------------------------------------------

if(! $opt_f ){
   print "Usage: $0\n";
   print "-f  sequences to scaffold (Multi-FASTA format, required)\n"; 
   print "-k  k-mer value (default -k $k, optional)\n";
   die "-p  Bloom filter false positive rate (default -p $fpr, optional - increase to prevent memory allocation errors)\n";
}

my $assemblyfile = $opt_f;
$k = $opt_k if($opt_k);
$fpr = $opt_p if($opt_p); 

print "\nRunning:$0 -f $assemblyfile -k $k -p $fpr\n\n";

my $bfout = $assemblyfile . "_k" . $k . "_p" . $fpr . "_rolling.bf";

if(! -e $assemblyfile){
   my $file_message = "\nInvalid file: $assemblyfile -- fatal\n";
   print $file_message;
   exit;
}else{
   my $file_message = "Checking sequence target file $assemblyfile...ok\n";
   print $file_message;
}


my $date = `date`;
chomp($date);
my $bloom;

eval{


print "$date:Estimating number of elements from file size\n";

my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$atime,$mtime,$ctime,$blksize,$blocks) = stat($assemblyfile);
my $bfelements = int($size);

my $date = `date`;
chomp($date);

#$bfelements *= 2;###twice the kmer when looking at both strands
my $m = ceil((-1 * $bfelements * log($fpr)) / (log(2) * log(2))); #the number of bits
my $rem = 64 - ($m % 64);
$m = $m + $rem;
my $hashfct = floor(($m / $bfelements) * log(2));# the number of hash functions
print "*****\nBloom filter specs\nelements=$bfelements\nFPR=$fpr\nsize (bits)=$m\nhash functions=$hashfct\n*****\n";
#$bloom = new Bloom::Faster({m=>$m,k=>$hashfct});
$bloom = new BloomFilter::BloomFilter($m, $hashfct, $k);

print "$date:Shredding supplied sequence file (-f $assemblyfile) into $k-mers..\n";
$bloom = &contigsToBloom($assemblyfile,$hashfct,$k,$bloom);

my $date = `date`;
chomp($date);


my $writing_tigbloom_message = "\n$date:Writing Bloom filter to disk ($bfout)\n";
print $writing_tigbloom_message;

$bloom->storeFilter($bfout);

#$bloom->to_file($bfout);###write BF to disk

};

if($@){
   my $message = $@;
   my $failure = "\nSomething went wrong running $0 $date\n$message\n";
   print $failure;
}else{
   my $success = "\n$date:$0 executed normally\n";
   print $success;
}

# manually free the data structures 
$bloom->DESTROY;

exit 1;


#----------------
sub contigsToBloom{
   my ($file,$hashfct,$k,$bloom) = @_;

   my $prevhead = "";
   my $seq = "";
   my $cttig=0;
   open(IN,$file) || die "Error reading $file -- fatal.\n";

   print "Contigs processed k=$k:\n";
   ###
   while(<IN>){
       chomp;
      if(/^\>(\S+)/){
         my $head=$1;

         if ($head ne $prevhead && $seq ne '' && $prevhead ne ''){
            $cttig++;
            print "\r$cttig";
            $|++;
            $bloom = &kmerizeContigBloom_newloop(uc($seq),$bloom,$hashfct,$k);
            #$bloom = &kmerizeContigBloom(uc($seq),$bloom,$k,$cttig,0);
            #my $revcomp = &reverseComplement(uc($seq));
            #$bloom = &kmerizeContigBloom($revcomp,$bloom,$k,$cttig,1);
         }
         $seq = '';
         $prevhead = $head;
      }else{
         $seq .= $_;
      }
   }
   $cttig++;
   print "\r$cttig";
   $|++;
   $bloom = &kmerizeContigBloom_newloop(uc($seq),$bloom,$hashfct,$k);
   #$bloom = &kmerizeContigBloom(uc($seq),$bloom,$k,$cttig,0);
   #my $revcomp = &reverseComplement(uc($seq));
   #$bloom = &kmerizeContigBloom($revcomp,$bloom,$k,$cttig,1);
   ###
   close IN;
   
   return $bloom;
}

#----------------
sub kmerizeContigBloom{
   my ($seq,$bloom,$hashfuct, $k) = @_;

   for(my $pos=0;$pos<=(length($seq)-$k);$pos++){
      my $kmer = substr($seq,$pos,$k);
      $bloom->insert($kmer);
   }
   return $bloom;
}

#sub kmerizeContigBloom_new{
#    my ($seq,$bloom,$hashfct,$k) = @_;
#
#    my $itr = new BloomFilter::ntHashIterator($seq, $hashfct, $k);
#    while($itr.neq(BloomFilter::ntHashIterator::end())) {
#        $bloom->insert($itr->deref());
#        $itr->incr();
#    }
#    return $bloom;
#}

sub kmerizeContigBloom_newloop{
    my ($seq,$bloom,$hashfct,$k) = @_;
    
    BloomFilter::insertSeq($bloom, $seq, $hashfct, $k);
    return $bloom
}

#-----------------------
sub reverseComplement{
   $_ = shift;
   $_ = uc();
   tr/ATGCYRKMBDHV/TACGRYMKVHDB/; 
   return (reverse());
}
            
