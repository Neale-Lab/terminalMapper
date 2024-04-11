#!/usr/bin/env perl
#Package Version: 2.0

######################################################################################################
# Author(s): T.J.Cooper
# Updated: 15/2/2017
# Trims FASTQ-entries from the 5'-end using user-specified patterns
######################################################################################################

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(basename);
use List::Util qw(all);
my $scriptname = basename($0);	#Obtain script-name
my ($R1,$R2,@UTrim,@LTrim);
my $usage = "Usage: $scriptname -1 <Read-1 FASTQ> -2 <Read-2 FASTQ> -U <5' Upstream Pattern> -L <3' Downstream Pattern>";
GetOptions('1=s' => \$R1,
           '2=s' => \$R2,
           'U=s' => \@UTrim,
           'L=s' => \@LTrim) or die("\n$usage\n");
die("\nError: Arguments or -flags are missing and/or incorrectly specified.\n\n$usage\n\n") unless all {defined} $R1, $R2, @UTrim, @LTrim;
my $outfile = "Trim_".basename($R1);
open my $OUT, '>', "$outfile" or die "$!";
open my $IN, '<', $R1 or die "$!";
my ($trimseq,$R1trimU,$R1trimL,$match,$R2trim,$platform);
while (<$IN>) {
   chomp $_;
   if ($. % 4 == 2) {
      $trimseq = $_;
      $platform = substr($trimseq,0,9);
      my $Uregex = join "|", @UTrim;
      if ($platform =~ m/^.*?$Uregex/) {
         $R1trimU = $+[0];
         $trimseq = substr($_,$R1trimU);
      } else {
         $R1trimU = 0;
      }
      my $Lregex = join "|", @LTrim;
      my $revseq = reverse($trimseq);
      if ($revseq =~ m/^.*?$Lregex/) {
         $R1trimL = $+[0];
         $revseq = substr($revseq,$R1trimL);
      } else {
         $revseq = substr($revseq,2);
         $R1trimL = 2;
      }
      my $finalseq = reverse($revseq);
      print $OUT "$finalseq\n";
   } elsif ($. % 4 == 0) {
      my $revqual = reverse($_);
      my $trimqual = substr($revqual,$R1trimL);
      $trimqual = reverse($trimqual);
      $trimqual = substr($trimqual,$R1trimU);
      print $OUT "$trimqual\n";
   } else {
      print $OUT "$_\n";
   }
}
my $outfile2 = "Trim_".basename($R2);
open my $OUT2, '>', "$outfile2" or die "$!";
open my $IN2, '<', $R2 or die "$!";
while (<$IN2>) {
   chomp $_;
   if ($. % 4 == 2) {
      $trimseq = $_;
      $platform = substr($trimseq,0,4);
      my $Uregex = join "|", @UTrim;
      if ($platform =~ m/^.*?$Uregex/) {
         $R1trimU = $+[0];
         $trimseq = substr($_,$R1trimU);
      } else {
         $R1trimU = 0;
      }
      my $Lregex = join "|", @LTrim;
      my $revseq = reverse($trimseq);
      if ($revseq =~ m/^.*?$Lregex/) {
         $R1trimL = $+[0];
         $revseq = substr($revseq,$R1trimL);
      } else {
         $revseq = substr($revseq,2);
         $R1trimL = 2;
      }
      my $finalseq = reverse($revseq);
      print $OUT2 "$finalseq\n";
   } elsif ($. % 4 == 0) {
      my $revqual = reverse($_);
      my $trimqual = substr($revqual,$R1trimL);
      $trimqual = reverse($trimqual);
      $trimqual = substr($trimqual,$R1trimU);
      print $OUT2 "$trimqual\n";
   } else {
      print $OUT2 "$_\n";
   }
}
