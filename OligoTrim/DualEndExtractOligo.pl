#!/usr/bin/env perl
#Package Version: 2.0

##################################################################################################################
# Author(s): T.J.Cooper
# Updated: 15/2/2017
# Processes paired-end .SAM files, extracting Watson + Crick coordinate information for Spo11 oligo libraries
# Calculates inter-event distances (between double-cut DSBs) and tallies instances of specific double-cuts
##################################################################################################################

use strict;
use warnings;
use Cwd;
use List::Util qw(first);
use Sort::Naturally;
local $|=1;
my $outext = '.txt';	#Output .file-extension
my $inext = $ARGV[0];	#Input .file-extension
my @files = glob("*$inext");
my $chk = scalar(@files);
print "\nFailed to detect any .SAM   files within the current directory.\n\n" if $chk == 0;
exit if $chk == 0;	#Stop script if no .SAM files are found
my $sub = cwd()."/Coordinates";
mkdir("$sub") unless $chk == 0;
my $sub2 = cwd()."/Analysis";
mkdir("$sub2") unless $chk == 0;
my $genome = $ARGV[2];
my $trimmode = $ARGV[3];
my $trimlength = $ARGV[4];
print "-------------------------------------";
print "\nCalculating Coordinates....\n";
print "-------------------------------------\n";
print "Currently processing:\n";
for my $file (@files) {   #For-each input file
	open my $IN, '<', $file or die "$!";	#Open and read input .SAM file(s)
	(my $strain = $file) =~ s/_[^_]+$//;	#Strain-name
	(my $mode = $ARGV[0]) =~ s/\.SAM//;		#Alignment-mode
	print "$strain\n";
	my $outfile = "Coordinates.".$genome."_".$strain."_".$mode.$outext;	#Output files
	my $outfile2 = "Coordinates.".$genome."_".$strain."_Ambiguous".$outext;
	open my $OUT, '>', "$sub/$outfile" or die "$!";
	open my $OUT2, '>', "$sub/$outfile2" or die "$!";
	print $OUT "PairID\tStrand\tChr\tPos\tReadLength\tCIGAR\tAdjustment\tRead-Flag\tmSize\n";
	print $OUT2 "PairID\tStrand\tChr\tPos\tReadLength\tCIGAR\tAdjustment\tRead-Flag\tmSize\tMD-Tag\t\n";
	my ($ID,$AID,$A,$B,$R1var,$R2var,$Wflag,$Cflag,%IED,%DoubleCut);
	while (<$IN>) {	#For-each .SAM record
		chomp $_;
		next if /^\s*@/; #Skip .SAM headerlines
		my @F = split("\t", $_);	#Split each tab-delimited field
		sub parseSAM {		#Subroutine to interpret SAM-field data
			my @rcd = @_;
			my @read;
			my $index = first{/MD:Z/} @rcd;	#Obtain variable-column MD:Z tag
			my @MDtag = $index =~ /\d+/g;		#Remove non-numeric characters
			my %rules = (M => 1,D => 1,I => 0,S => 1);	#Rules to handle insertion/deletions/matches/soft-clipping
			my ($s,$LS,$RS) = (0)x3;
			while ($rcd[5] =~ /(\d+)([MDIS])/g) {		#Parse and interpret CIGAR code
				my ($n,$op)  = ($1,$2);
				$s += $n * $rules{$op} unless $op eq 'S';		#Calculate POS adjustment (insertions/deletions)
				$LS += $n * $rules{$op} if $op eq 'S' && $-[0]==0; 	#(upstream soft-clip)
				$RS += $n * $rules{$op} if $op eq 'S' && $+[0]==length($rcd[5]); 	#(downstream soft-clip)
			}
			my $l = length($rcd[9]);	#Read-length
			my $wp = $rcd[3]-$LS;	 #Adjusted 5' coordinate (Watson strand)
			my $cp = $rcd[3]+($RS+$s)-1;	 #Adjusted 5' coordinate (Crick strand)
			push(@read, $rcd[2],$wp,$cp,$l,$rcd[5],$s,$LS,$RS,$index);
			return(\@read, \@MDtag);
		}
		if (grep {$_ == $F[1]} 99,83) {	#For 99/147 or 83/163 read-pairs
			my $partner = <$IN>;
			my @F2 = split("\t", $partner);	#Split each tab-delimited field
			if ($F[1] == 99) {
				($A, $R1var) = parseSAM(@F); ($B, $R2var) = parseSAM(@F2);
			} else {
				($A, $R1var) = parseSAM(@F2); ($B, $R2var) = parseSAM(@F);
			}
			my @revMDtag = reverse(@{$R2var});
			my @Wat = @{$A}; my @Cri = @{$B};
			my ($chr, $pos, $rl, $cigar, $Lclip, $vtag) = @Wat[0,1,3,4,6,8];
			my ($chrp, $posp, $rlp, $cigarp, $cc, $Rclip, $vtagp) = @Cri[0,2,3,4,5,7,8];
			if ($F[1] == 99) {
				$Wflag = 99; $Cflag = 147; $posp = $posp+2;
				$DoubleCut{$chr}{$pos}{$posp}{"99"}++;
			} elsif ($F[1] == 83) {
				$Wflag = 163; $Cflag = 83; $pos = $pos-2;
				$DoubleCut{$chr}{$pos}{$posp}{"163"}++;
			}
			my $mlength = abs($pos-$posp)+1;
			$IED{$mlength}++;	#Calculates and tallies IEDs
			if ($R1var->[0] == 0 && $R1var->[1] == 0 && $pos > 0 || $revMDtag[0] == 0 && $revMDtag[1] == 0 && $pos > 0 || $Lclip > 1 && $pos > 0 || $Rclip > 1 && $pos > 0) {	#Detect ambigious ends (99/147 or 83/163 pairs)
				$AID++;
				printf($OUT2 "%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%s\n%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%s\n", $AID,"w",$chr,$pos,$rl,$cigar,0-$Lclip,$Wflag,$mlength,$vtag,$AID,"c",$chrp,$posp,$rlp,$cigarp,$Rclip+$cc-1,$Cflag,$mlength,$vtagp);
			} else {
				$ID++;	#Detect ambigious ends (99/147 or 83/163 pairs)
				printf($OUT "%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\n%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\n", $ID,"w",$chr,$pos,$rl,$cigar,0-$Lclip,$Wflag,$mlength,$ID,"c",$chrp,$posp,$rlp,$cigarp,$Rclip+$cc-1,$Cflag,$mlength);
			}
		}
	}
	close $IN;
	close $OUT;
	my $outfile3 = "IED.".$genome."_".$strain.$outext;
	my $outfile4 = "DoubleCuts.".$genome."_".$strain.$outext;
	open my $OUT3, '>', "$sub2/$outfile3" or die "$!";
	open my $OUT4, '>', "$sub2/$outfile4" or die "$!";
	print $OUT3 "IED\tFreq\n";
	print $OUT4 "Chr\tCoord-A\tCoord-B\tR1W Freq\tR1C Freq\n";
	foreach my $key (sort {$a <=> $b} keys %IED) {
		print $OUT3 "$key\t$IED{$key}\n";
	}
	foreach my $chrn (nsort keys %DoubleCut) {
		foreach my $coord1 (sort {$a <=> $b} keys %{$DoubleCut{$chrn}}) {
			foreach my $coord2 (sort {$a <=> $b} keys %{$DoubleCut{$chrn}{$coord1}}) {
					print $OUT4 join("\t",$chrn, $coord1,$coord2,map $_ // 0, @{$DoubleCut{$chrn}{$coord1}{$coord2}}{qw{99 163}});
					print $OUT4 "\n";
			}
		}
	}
	if ($ARGV[1] eq "Y") {
		chomp($file);
		unlink($file);
	}
}
