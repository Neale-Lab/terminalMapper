#!/usr/bin/env perl
#Package Version: 3.0

##############################################################################################################################
# Author(s): T.J.Cooper
# Updated: 6/12/2016
# Generates 1bp-histogram (sparse format) from combined coordinate files
##############################################################################################################################

use strict;
use warnings;
use Cwd;
use Getopt::Long;
use File::Basename qw(basename);
use List::Util qw(all);
use Sort::Naturally;
my $scriptname = basename($0);	#Obtain script-name
my ($coord,$genome,$target);
my $usage = "Usage: $scriptname -i <Input Folder> -g <Genome Name>";
GetOptions('i=s' => \$coord,
           'g=s' => \$genome) or die("\n$usage\n");
die("\nError: Arguments or -flags are missing and/or incorrectly specified.\n\n$usage\n\n") unless all {defined} $coord, $genome;
my @files = glob($coord."/*_Combined.txt");
my $chk = scalar(@files);
if ($chk == 0) {
   @files = glob($coord."/*_Global.txt");
   $chk = scalar(@files);
}
$genome =~ s/_//g;
print "\nFailed to detect any valid coordinate files within the specified directory.\n\n" if $chk == 0;
exit if $chk == 0;
for my $file (@files) {
	open my $IN, '<', $file or die "$!";
	(my $strain = basename($file)) =~ s/_[^_]+$//;	#Strain-name
   $strain =~ s/^[^_]*_//;
	(my $wd = $coord) =~ s/[^\/]+$//;
	my (%hits,%R1hits,%R2hits);
	my @headers = split("\t",<$IN>);
	my %headidx; @headidx{@headers} = 0 .. $#headers;
	my $index = $headidx{Strand};
	while (<$IN>) {
		chomp $_;
		my @F = split("\t", $_);
		$hits{$F[$index+1]}{$F[$index+2]}{$F[$index]}++;
      if ($F[$index+6] == 99 || $F[$index+6] == 83) {
         $R1hits{$F[$index+1]}{$F[$index+2]}{$F[$index]}++;
      }
      if ($F[$index+6] == 147 || $F[$index+6] == 163) {
         $R2hits{$F[$index+1]}{$F[$index+2]}{$F[$index]}++;
      }
	}
   my $outfile = $wd."/Histograms/"."Map.R1R2_".$genome."_".$strain.".txt";
   open my $OUT, '>', "$outfile" or die "$!";
   print $OUT "Chr\tPos\tWatson\tCrick\n";
	foreach my $chr (nsort keys %hits) {
	    foreach my $pos (sort {$a <=> $b} keys %{ $hits{$chr} }) {
	        print $OUT join("\t", $chr, $pos, map $_ // 0, @{ $hits{$chr}{$pos} }{qw{w c}});
			  print $OUT "\n";
	    }
	}
   if (%R2hits) {
      my $outfile2 = $wd."/Histograms/"."Map.R1_".$genome."_".$strain.".txt";
      my $outfile3 = $wd."/Histograms/"."Map.R2_".$genome."_".$strain.".txt";
      open my $OUT2, '>', "$outfile2" or die "$!";
      open my $OUT3, '>', "$outfile3" or die "$!";
      for my $OF ($OUT2, $OUT3) { print $OF "Chr\tPos\tWatson\tCrick\n"; }
      foreach my $chr (nsort keys %R1hits) {
          foreach my $pos (sort {$a <=> $b} keys %{ $R1hits{$chr} }) {
              print $OUT2 join("\t", $chr, $pos, map $_ // 0, @{ $R1hits{$chr}{$pos} }{qw{w c}});
              print $OUT2 "\n";
          }
      }
      foreach my $chr (nsort keys %R2hits) {
          foreach my $pos (sort {$a <=> $b} keys %{ $R2hits{$chr} }) {
              print $OUT3 join("\t", $chr, $pos, map $_ // 0, @{ $R2hits{$chr}{$pos} }{qw{w c}});
              print $OUT3 "\n";
          }
      }
   }
}
