#!/usr/bin/env perl
#Package Version: 4.0

##################################################################################################################
# Author(s): T.J.Cooper
# Updated: 1/2/2018
# Processes paired-end .SAM files, extracting Watson + Crick coordinate information for double-cut Spo11 libraries
# Quality-control and filtering (atypical read-orientation, dubious ends)
# Two-step alignment (unmapped mate read-trimming, --local alignment)
# Calculates molecule sizes and tallies instances of specific double-cuts
##################################################################################################################

use strict;
use warnings;
use Cwd;
use List::Util qw(first any);
use Sort::Naturally;
use File::Basename;
local $|=1;
my ($name, $path, $suffix) = fileparse($0);	#Obtain script-name and directory
my $outExt = '.txt';	#Output .file-extension
my $alignRound = $ARGV[0];	#Input .file-extension
my @files = glob("*$alignRound");
my $fileCheck = scalar(@files);
print "\nFailed to detect any .SAM files within the current directory.\n\n" if $fileCheck == 0;
exit if $fileCheck == 0;	#Stop script if no .SAM files are found
my $subDir = cwd()."/Data";
mkdir("$subDir");
my $subDir2 = cwd()."/Analysis";
mkdir("$subDir2");
my $genomeName = $ARGV[2];
my $trimMode = $ARGV[3];
my $trimLength = $ARGV[4];
my (%list,$chr,$min,$max,%chrlabels,%fh);
open my $IN, '<', $ARGV[8] or die "$!";
while (<$IN>) {
	chomp $_;
	my @F = split("\t", $_);
	$chrlabels{$F[0]} = undef;
}
if (defined $ARGV[9]) {
	open my $IN2, '<', $ARGV[9] or die "$!";
	<$IN2> for (1..1);
	while (<$IN2>) {
		chomp $_;
		($chr, $min, $max) = split('\s+', $_);
		push @{$list{$chr}}, [$min,$max];
	}
}
print "-------------------------------------";
print "\nCalculating Coordinates....\n";
print "-------------------------------------\n";
print "Currently processing:\n";
my $outfile5 = cwd()."/Logs/RunSummary".$outExt;
open my $OUT5, '>>', "$outfile5" or die "$!";
print $OUT5 "ALIGNMENT & MAPPING","\t" x (8),"CALLING","\n";
print $OUT5 join("\t","Sample","Total Read Pairs (A)","Total Mapped Pairs (B)","% of (A)","Mapped Pairs (Unique)","% of (B)","Mapped Pairs (Multi)","% of (B)","Valid Pairs (C)","% of (B)","Ambiguous Pairs","% of (B)","Repeat Pairs","% of (C)"),"\n";
for my $file (@files) {   #For-each input file
	open my $IN3, '<', $file or die "$!";	#Open and read input .SAM file(s)
	(my $strain = $file) =~ s/_[^_]+$//;	#Strain-name
	(my $mode = $ARGV[0]) =~ s/\.SAM//;		#Alignment-mode
	mkdir("$subDir/$strain");
	mkdir("$subDir2/$strain");
	print "$strain\n";
	for my $label (keys %chrlabels) {
		my $outfile = "Data.".$label."_".$genomeName."_".$strain."_".$mode.$outExt;	#Output files
		open $fh{$label}, '>', "$subDir/$strain/$outfile" or die "$!";
		print {$fh{$label}} "PairID\tStrand\tChr\tPos\tReadLength\tCIGAR\tAdjustment\tRead-Flag\tmSize\tRepeat\n";
	}
	my $outfile2 = "Data.".$genomeName."_".$strain."_Ambiguous".$outExt;
	my ($OUT, $OUT2, $OUT3, $OUT4);
	open $OUT2, '>>', "$subDir/$strain/$outfile2" or die "$!";
	print $OUT2 "PairID\tStrand\tChr\tPos\tReadLength\tCIGAR\tAdjustment\tRead-Flag\tmSize\tMD-Tag\n";
	if ($alignRound eq "Global.SAM" && $trimMode eq "Y") {
		my $outfile3 = $strain.$ARGV[5]."_unmapped_trimmed.fastq";	#Unmapped R1 FASTQ file
		my $outfile4 = $strain.$ARGV[6]."_unmapped_trimmed.fastq";	#Unmapped R2 FASTQ file
		open $OUT3, '>', "$outfile3" or die "$!";
		open $OUT4, '>', "$outfile4" or die "$!";
	}
	my ($ID,$AID,$A,$B,$R1var,$R2var,$Wflag,$Cflag,%Size,%DoubleCut,$multimap);
	my $repeat = (0);
	while (<$IN3>) {	#For-each .SAM record
		chomp $_;
		next if /^\s*@/; #Skip .SAM headerlines
		my @F = split("\t", $_);	#Split each tab-delimited field
		my $orientation = $F[3]-$F[7];	#Discard atypical read-orientations
		if ($F[1] == 99 && $orientation > 0 || $F[1] == 83 && $orientation < 0) {
			my $skipline = <$IN3>;
			next;
		}
		unless (exists($chrlabels{$F[2]})) {
			my $skipline = <$IN3>;
			next;
		}
		if ($alignRound eq "Global.SAM" && $trimMode eq "Y") {	#Populate unmapped R1/R2 FASTQ files mapped-unmapped pairs
			if (grep {$_ == $F[1]} 73,137) {
				print $OUT3 "\@$F[0] 1:N:0:1\n$F[9]\n+\n$F[10]\n" if $F[1] == 73;
				print $OUT4 "\@$F[0] 1:N:0:1\n$F[9]\n+\n$F[10]\n" if $F[1] == 137;
			}
			if (grep {$_ == $F[1]} 89,153) {
				$F[9] =~ tr/GATC/CTAG/;
				my $revseq = reverse($F[9]);
				my $revqual = reverse($F[10]);
				print $OUT3 "\@$F[0] 1:N:0:1\n$revseq\n+\n$revqual\n" if $F[1] == 89;
				print $OUT4 "\@$F[0] 1:N:0:1\n$revseq\n+\n$revqual\n" if $F[1] == 153;
			}
			if (grep {$_ == $F[1]} 69,133) {
				my $trimseq = subDirstr($F[9],0,$trimLength);
				my $trimqual = subDirstr($F[10],0,$trimLength);
				print $OUT3 "\@$F[0] 1:N:0:1\n$trimseq\n+\n$trimqual\n" if $F[1] == 69;
				print $OUT4 "\@$F[0] 1:N:0:1\n$trimseq\n+\n$trimqual\n" if $F[1] == 133;
			}
		}
		sub parseSAM {		#subDirroutine to interpret SAM-field data
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
		my $rflag = 0;
		if (grep {$_ == $F[1]} 99,83) {	#For 99/147 or 83/163 read-pairs
			my $partner = <$IN3>;
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
			my $mlength = abs($pos-$posp);
			if ($F[1] == 99) {
				$Wflag = 99; $Cflag = 147;
				if ($ARGV[0] eq "Global.SAM") {
					$Size{$mlength}++;
				}
			} elsif ($F[1] == 83) {
				$Wflag = 163; $Cflag = 83;
				if ($ARGV[0] eq "Global.SAM") {
					$Size{$mlength}++;
				}
			}
			if ($R1var->[0] == 0 && $R1var->[1] == 0 && $pos > 0 || $revMDtag[0] == 0 && $revMDtag[1] == 0 && $pos > 0 && $ARGV[7] eq "DOUBLE" || $Lclip > 1 && $pos > 0 || $Rclip > 1 && $pos > 0 && $ARGV[7] eq "DOUBLE") {	#Detect ambigious ends (99/147 or 83/163 pairs)
				$AID++;
				print $OUT2 join("\t",$AID,"w",$chr,$pos,$rl,$cigar,0-$Lclip,$Wflag,$mlength,$vtag),"\n",join("\t",$AID,"c",$chrp,$posp,$rlp,$cigarp,$Rclip+$cc-1,$Cflag,$mlength,$vtagp),"\n";
			} else {
				$ID++;
				if (grep {$pos >= $_->[0] && $pos <= $_->[1] || $posp >= $_->[0] && $posp <= $_->[1]} @{ $list{$chr}}) {
					$repeat++; $rflag = 1;
				}
				print {$fh{$F[2]}} join("\t",$ID,"w",$chr,$pos,$rl,$cigar,0-$Lclip,$Wflag,$mlength,$rflag),"\n",join("\t",$ID,"c",$chrp,$posp,$rlp,$cigarp,$Rclip+$cc-1,$Cflag,$mlength,$rflag),"\n";
			}
		}
	}
	close $IN3;
	close $OUT2;
	my $log = cwd()."/Logs/Log.".$strain.$outExt;
	my ($TRP,$MP,$MM,$TM,$AR,$MMR,$CPR,$ARR);
	open my $IN4, '+<', $log or die "$!";
	while (<$IN4>) {	#For-each .SAM record
		chomp $_;
		if (/(\d+) \(\d+.\d+\%\) were paired/) { $TRP = $1; }
		if (/(\d+) \(\d+.\d+\%\) aligned concordantly exactly/) { $MP = $1; }
		if (/(\d+) \(\d+.\d+\%\) aligned concordantly >/) { $MM = $1; }
	}
	print $IN4 "\n---------------------------------------\nCALL STATS\n---------------------------------------\n";
	print $IN4 "Total Pairs:\t",$TRP,"\n";
	printf($IN4 "%s\t%d\n","Valid Pairs:",$ID);
	printf($IN4 "%s\t%d\n","Repeat Pairs:",$repeat);
	printf($IN4 "%s\t%d\n","Ambig Pairs:",$AID);
	if ($ARGV[0] eq "Global.SAM") {
		my $totalMapped = $MM+$MP;
		my $unMapped = $TRP-($totalMapped);
		my $uniqueAdjust = $MP-$MM;
		my $validAdjust = $ID-($repeat+$AID);
		my $repeatAdjust = $repeat-$AID;
		$AR = (int(($MP/$TRP)*10000)/100); $MMR = (int(($MM/$TRP)*10000)/100);
		$ARR = (int(($AID/($MP+$MM))*10000)/100); my $DR = (int(($repeat/($MP+$MM))*10000)/100);
		$CPR = (int(($ID/($MP+$MM))*10000)/100)-$DR;
		my $uniqueRate = (int(($MP/($MP+$MM))*10000)/100);
		my $multiRate = (int(($MM/($MP+$MM))*10000)/100);
		print $OUT5 join("\t",$strain,$TRP,$MP+$MM,$MMR+$AR,$MP,$uniqueRate,$MM,$multiRate,$ID,$CPR,$AID,$ARR,$repeat,$DR),"\n";
		my $outfile6 = "/".$strain."/"."Sizes.".$genomeName."_".$strain.$outExt;
		my $outfile7 = "Coordinates.".$genomeName."_".$strain.$outExt;
		open my $OUT6, '>', "$subDir2/$outfile6" or die "$!";
		print $OUT6 "Size\tFreq\n";
		foreach my $key (sort {$a <=> $b} keys %Size) {
			print $OUT6 "$key\t$Size{$key}\n";
		}
		my $Routfile = "Plot.".$strain;
		my $RScript = $path."alignPlots.R";
		my $command = "Rscript '$RScript' $strain $unMapped $totalMapped $MP $MM $validAdjust $repeat $AID $Routfile";
		system($command);
	}
	if ($ARGV[1] eq "Y") {
		chomp($file);
		unlink($file);
	}
}
