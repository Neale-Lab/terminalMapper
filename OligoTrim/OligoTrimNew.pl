#!/usr/bin/env perl
#Package Version: 2.0

######################################################################################################
# Author(s): T.J.Cooper
# Updated: 19/6/2018
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
my $outfile2 = "Trim_".basename($R2);
open my $OUT, '>', "$outfile" or die "$!";
open my $IN, '<', $R1 or die "$!";
open my $OUT2, '>', "$outfile2" or die "$!";
open my $IN2, '<', $R2 or die "$!";
my ($trimseq_R1,$trimseq_R2,$R1trimU,$R1trimL,$match,$R2trimU,$R2trimL,$platform_R1,$platform_R2);
while (<$IN>) {
  chomp;
  my $headerline_R1 = $_; my $headerline_R2 = <$IN2>;
  my $seqline_R1 = <$IN>; my $seqline_R2 = <$IN2>;
  my $infoline_R1 = <$IN>; my $infoline_R2 = <$IN2>;
  my $qline_R1 = <$IN>; my $qline_R2 = <$IN2>;
  $platform_R1 = substr($seqline_R1,0,9);
  $platform_R2 = substr($seqline_R2,0,4);
  my $checkR1 = 0; my $checkR2 = 0;
  foreach my $A1 (@UTrim) {
    if ($platform_R1 =~ m/^.*?$A1/) {
       $checkR1 = 1;
       $R1trimU = $+[0];
       $trimseq_R1 = substr($seqline_R1,$R1trimU);
       last;
     }
  }
  foreach my $A2 (@UTrim) {
    if ($platform_R2 =~ m/^.*?$A2/) {
       $checkR2 = 1;
       $R2trimU = $+[0];
       $trimseq_R2 = substr($seqline_R2,$R2trimU);
       last;
     }
  }
  if ($checkR1==1 and $checkR2==1) {
    $R1trimL = 0;
    foreach my $B1 (@LTrim) {
      if ($trimseq_R1 =~ m/^.*?$B1/) {
        $R1trimL = length($trimseq_R1)-$+[0]+length($B1);
        my $revseq_R1 = reverse($trimseq_R1);
        $trimseq_R1 = substr($revseq_R1,$R1trimL);
        last;
      }
    }
    $R2trimL = 0;
    foreach my $B2 (@LTrim) {
      if ($trimseq_R2 =~ m/^.*?$B2/) {
        $R2trimL = length($trimseq_R2)-$+[0]+length($B2);
        my $revseq_R2 = reverse($trimseq_R2);
        $trimseq_R2 = substr($revseq_R2,$R2trimL);
        last;
      }
    }
    $trimseq_R1 = reverse($trimseq_R1) unless $R1trimL==0;
    my $revqual_R1 = reverse($qline_R1);
    my $trimqual_R1 = substr($revqual_R1,$R1trimL);
    $trimqual_R1 = reverse($trimqual_R1);
    $trimqual_R1 = substr($trimqual_R1,$R1trimU);
    print $OUT "$headerline_R1\n$trimseq_R1\n$infoline_R1$trimqual_R1\n";
    $trimseq_R2 = reverse($trimseq_R2) unless $R2trimL==0;
    my $revqual_R2 = reverse($qline_R2);
    my $trimqual_R2 = substr($revqual_R2,$R2trimL);
    $trimqual_R2 = reverse($trimqual_R2);
    $trimqual_R2 = substr($trimqual_R2,$R2trimU);
    print $OUT2 "$headerline_R2$trimseq_R2\n$infoline_R2$trimqual_R2\n";
  }
}
