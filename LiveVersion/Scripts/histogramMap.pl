#!/usr/bin/env perl
#Package Version: 4.1

##############################################################################################################################
# Author(s): T.J.Cooper
# Updated: 1/2/2018
# Generates 1bp-histogram (sparse format) from combined coordinate files
##############################################################################################################################

use strict;
use warnings;
use Cwd;
use Getopt::Long;
use File::Basename qw(basename);
use List::Util qw(all);
use Sort::Naturally;
my $scriptname = basename($0);  #Obtain script-name
my ($coord,$genome,$target);
my $usage = "Usage: $scriptname -i <Input Folder> -g <Genome Name>";
GetOptions('i=s' => \$coord,
           'g=s' => \$genome) or die("\n$usage\n");
die("\nError: Arguments or -flags are missing and/or incorrectly specified.\n\n$usage\n\n") unless all {defined} $coord, $genome;
my @files = glob($coord."/*/*_Calls.txt");
my $chk = scalar(@files);
print "\nFailed to detect any valid coordinate files within the specified directory.\n\n" if $chk == 0;
exit if $chk == 0;
$genome =~ s/_//g;
for my $file (@files) {
	open my $IN, '<', $file or die "$!";
	(my $strain = basename($file)) =~ s/_Calls.txt//;	#Strain-name
  (my $prefix = $strain) =~ s/Data.//;
  $strain =~ s/[^_]*_[^_]*_//;
  mkdir("Histograms/$strain") unless (-d "Histograms/$strain");
	(my $wd = $coord) =~ s/[^\/]+$//;
	my (%hits,%R1hits,%R2hits,%Duplicates);
	my @headers = split("\t",<$IN>);
	my %headidx; @headidx{@headers} = 0 .. $#headers;
	my $index = $headidx{Strand};
	while (<$IN>) {
		chomp $_;
		my @F = split("\t", $_);
    my $partner = <$IN>;
    my @F2 = split("\t", $partner);
    $hits{$F[$index+1]}{$F[$index+2]}{$F[$index]}++;  #Tally R1+R2 hits
    $hits{$F2[$index+1]}{$F2[$index+2]}{$F2[$index]}++;
    if ($F[$index+6] == 99) {
       $R1hits{$F[$index+1]}{$F[$index+2]}{$F[$index]}++; #Tally R1 (99/83) hits
       $R2hits{$F2[$index+1]}{$F2[$index+2]}{$F2[$index]}++;  #Tally R2 (147/163) hits
       $Duplicates{$F[2]}{$F[3]}{$F2[3]}{"99"}++; #Tally 5'-3' molecule coordinates
    }
    if ($F[$index+6] == 163) {
       $R1hits{$F2[$index+1]}{$F2[$index+2]}{$F2[$index]}++;
       $R2hits{$F[$index+1]}{$F[$index+2]}{$F[$index]}++;
       $Duplicates{$F[2]}{$F[3]}{$F2[3]}{"163"}++;
    }
	}
  my $outfile = $wd."Histograms/".$strain."/Map.R1R2_".$prefix.".txt";  #Print sparsely formatted histogram (total)
  open my $OUT, '>', "$outfile" or die "$!";
  print $OUT "Chr\tPos\tWatson\tCrick\n";
	foreach my $chr (nsort keys %hits) {
	    foreach my $pos (sort {$a <=> $b} keys %{ $hits{$chr} }) {
	        print $OUT join("\t", $chr, $pos, map $_ // 0, @{ $hits{$chr}{$pos} }{qw{w c}});
			  print $OUT "\n";
	    }
	}
  my $outfile2 = $wd."/Histograms/".$strain."/Map.R1_".$prefix.".txt";
  my $outfile3 = $wd."Histograms/".$strain."/Map.R2_".$prefix.".txt";
  open my $OUT2, '>', "$outfile2" or die "$!";
  open my $OUT3, '>', "$outfile3" or die "$!";
  for my $OF ($OUT2, $OUT3) { print $OF "Chr\tPos\tWatson\tCrick\n"; }
    foreach my $chr (nsort keys %R1hits) {  #Print sparsely formatted histogram (R1)
        foreach my $pos (sort {$a <=> $b} keys %{ $R1hits{$chr} }) {
            print $OUT2 join("\t", $chr, $pos, map $_ // 0, @{ $R1hits{$chr}{$pos} }{qw{w c}});
            print $OUT2 "\n";
        }
    }
    foreach my $chr (nsort keys %R2hits) {  #Print sparsely formatted histogram (R2)
        foreach my $pos (sort {$a <=> $b} keys %{ $R2hits{$chr} }) {
            print $OUT3 join("\t", $chr, $pos, map $_ // 0, @{ $R2hits{$chr}{$pos} }{qw{w c}});
            print $OUT3 "\n";
        }
    }

   my $outfile4 = $wd."Analysis/".$strain."/Molecules.".$prefix.".txt";
   open my $OUT4, '>', "$outfile4" or die "$!";
   print $OUT4 "Chr\tCoord-A\tCoord-B\tR1W Freq\tR1C Freq\n";
   foreach my $chrn (nsort keys %Duplicates) {  #Print sparsely formatted histogram (molcule sizes)
     foreach my $coord1 (sort {$a <=> $b} keys %{$Duplicates{$chrn}}) {
       foreach my $coord2 (sort {$a <=> $b} keys %{$Duplicates{$chrn}{$coord1}}) {
           print $OUT4 join("\t",$chrn,$coord1,$coord2,map $_ // 0, @{$Duplicates{$chrn}{$coord1}{$coord2}}{qw{99 163}});
           print $OUT4 "\n";
       }
     }
   }
}
