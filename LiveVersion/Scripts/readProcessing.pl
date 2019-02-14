#!/usr/bin/env perl
#Package Version: 4.1

########################################################################################################
# Author(s): T.J.Cooper
# Updated: 17/5/2018
# Processes paired-end .SAM files, extracting Watson + Crick coordinate information for ccSeq-libraries
# Quality-control and filtering (atypical read-orientation, dubious ends)
# Calculates molecule sizes (tally)
########################################################################################################

use strict;
use warnings;
use Cwd;
use List::Util qw(first any);
use File::Basename;
local $|=1;	#Auto-flush (buffering)

###########################
# Input Files
###########################
my ($name, $path, $suffix) = fileparse($0);	#Obtain script-name and directory
my $outExt = '.txt';	#Output .file-extension
my $file = $ARGV[0];
my $subDir = cwd()."/Data";	#Create output subdirectories
my $subDir2 = cwd()."/Analysis";
mkdir("$subDir");
mkdir("$subDir2");
my $genomeName = $ARGV[1];
my (%list,$chr,$min,$max,%chrlabels,%fh);
open my $IN, '<', $ARGV[6] or die "$!";	#Read and parse FASTA index
while (<$IN>) {
	chomp $_;
	my @F = split("\t", $_);
	$chrlabels{$F[0]} = undef;
}
if (defined $ARGV[8]) {	#Read and parse repeat list
	open my $IN2, '<', $ARGV[8] or die "$!";
	<$IN2> for (1..1);
	while (<$IN2>) {
		chomp $_;
		($chr, $min, $max) = split('\s+', $_);
		push @{$list{$chr}}, [$min,$max];
	}
}

###########################
# Log File
###########################
my $outfile = cwd()."/Logs/RunSummary".$outExt;
open my $OUT, '>>', "$outfile" or die "$!";
print $OUT "ALIGNMENT & MAPPING","\t" x (8),"CALLING","\n";
print $OUT join("\t",
    "Sample",
    "Total Read Pairs (A)",
    "Total Mapped Pairs (B)","% of (A)","Mapped Pairs (Unique)","% of (B)",
    "Mapped Pairs (Multi)","% of (B)","Valid Pairs (C)","% of (B)",
    "Ambiguous Pairs","% of (B)","Repeat Pairs","% of (C)"),"\n";

###########################
# Output Files
###########################
open my $IN3, '<', $file or die "$!";	#Open and read input .SAM file
(my $strain = $file) =~ s/_[^_]+$//;	#Strain-name
$strain =~ s/.\///;
print "$strain\n";
mkdir("$subDir/$strain");
mkdir("$subDir2/$strain");
for my $label (keys %chrlabels) {
	my $outfile2 = "Data.".$label."_".$genomeName."_".$strain."_Calls".$outExt;	#Output files
	open $fh{$label}, '>', "$subDir/$strain/$outfile2" or die "$!";
	print {$fh{$label}} join("\t","PairID","Strand","Chr","Pos","ReadLength","CIGAR","Adjustment","Read-Flag","mSize","Repeat","\n");
}
my $outfile3 = "Data.".$genomeName."_".$strain."_Ambiguous".$outExt;
open my $OUT2, '>>', "$subDir/$strain/$outfile3" or die "$!";
print $OUT2 join("\t","PairID","Strand","Chr","Pos","ReadLength","CIGAR","Adjustment","Read-Flag","mSize","MD-Tag","\n");

###########################
# Read Processing
###########################

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

	my $rflag = 0;
	if (grep {$_ == $F[1]} 99,83) {	#For 99/147 or 83/163 read-pairs
		my $partner = <$IN3>;
		if ($F[4] < $ARGV[7]) {
			next;
		}
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
			$Wflag = 99; $Cflag = 147;
		} elsif ($F[1] == 83) {
			$Wflag = 163; $Cflag = 83;
		}
		my $mlength = abs($pos-$posp)+1;
		$Size{$mlength}++;

		#Detect ambigious ends (99/147 or 83/163 pairs)
		if ($R1var->[0] == 0 && $R1var->[1] == 0 && $pos > 0 || #5' SNPs
		$revMDtag[0] == 0 && $revMDtag[1] == 0 && $pos > 0 && $ARGV[6] eq "DOUBLE" #3' SNPs (Double-cut library ONLY)
		|| $Lclip > 1 && $pos > 0	#5' --local clipping
		|| $Rclip > 1 && $pos > 0 && $ARGV[6] eq "DOUBLE") {	#3' --local clipping (Double-cut library ONLY)
			$AID++;
			print $OUT2 join("\t",$AID,"w",$chr,$pos,$rl,$cigar,0-$Lclip,$Wflag,$mlength,$vtag),	#Print ambigious reads
			"\n",join("\t",$AID,"c",$chrp,$posp,$rlp,$cigarp,$Rclip+$cc-1,$Cflag,$mlength,$vtagp),"\n";
		} else {	#Non-ambigious ends
			$ID++;
			if (grep {$pos >= $_->[0] && $pos <= $_->[1] || $posp >= $_->[0] && $posp <= $_->[1]} @{ $list{$chr}}) { #Detect reads mapping to repeat regions
				$repeat++; $rflag = 1;
			}
			print {$fh{$F[2]}} join("\t",$ID,"w",$chr,$pos,$rl,$cigar,0-$Lclip,$Wflag,$mlength,$rflag),	#Print non-ambigious reads
			"\n",join("\t",$ID,"c",$chrp,$posp,$rlp,$cigarp,$Rclip+$cc-1,$Cflag,$mlength,$rflag),"\n";
		}
	}
}
close $IN3;
close $OUT2;

###########################
# Alignment/Library Stats
###########################
my $log = cwd()."/Logs/Log.".$strain.$outExt;
my ($TRP,$MP,$MM,$TM,$AR,$MMR,$CPR,$ARR);
open my $IN4, '+<', $log or die "$!";	#Open and read Bowtie2 .log file
while (<$IN4>) {	#For-each .log record
	chomp $_;
	if (/(\d+) \(\d+.\d+\%\) were paired/) { $TRP = $1; }
	if (/(\d+) \(\d+.\d+\%\) aligned concordantly exactly/) { $MP = $1; }
	if (/(\d+) \(\d+.\d+\%\) aligned concordantly >/) { $MM = $1; }
}
print $IN4 "\n---------------------------------------\nCALL STATS\n---------------------------------------\n";
print $IN4 "Total Pairs:\t",$TRP,"\n";
printf($IN4 "%s\t%d\n%s\t%d\n%s\t%d\n","Valid Pairs:",$ID,"Repeat Pairs:",$repeat,"Ambig Pairs",$AID);
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
print $OUT join("\t",$strain,$TRP,$MP+$MM,$MMR+$AR,$MP,$uniqueRate,$MM,$multiRate,$ID,$CPR,$AID,$ARR,$repeat,$DR),"\n";
my $outfile6 = "/".$strain."/"."Sizes.".$genomeName."_".$strain.$outExt;
open my $OUT3, '>', "$subDir2/$outfile6" or die "$!";
print $OUT3 "Size\tFreq\n";	#Print molecule/ICD sizes
foreach my $key (sort {$a <=> $b} keys %Size) {
	print $OUT3 "$key\t$Size{$key}\n";
}

###########################
# Alignment Plots
###########################
my $Routfile = "Plot.".$strain;	#R-Script
my $RScript = $path."alignPlots.R";
my $command = "Rscript '$RScript' $strain $unMapped $totalMapped $MP $MM $validAdjust $repeat $AID $Routfile";
system($command);

###########################
# SpaceSaver
###########################
if ($ARGV[2] eq "Y") {	#Delete .SAM file if SPACE_SAVER=Y
	chomp($file);
	unlink($file);
}
