#!/usr/bin/env perl
#Package Version: 3.0

##################################################################################################################
# Author(s): T.J.Cooper
# Updated: 15/2/2017
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
use Chart::Gnuplot;
local $|=1;
my $outExt = '.txt';	#Output .file-extension
my $alignRound = $ARGV[0];	#Input .file-extension
my @files = glob("*$alignRound");
my $fileCheck = scalar(@files);
print "\nFailed to detect any .SAM files within the current directory.\n\n" if $fileCheck == 0;
exit if $fileCheck == 0;	#Stop script if no .SAM files are found
my $subDir = cwd()."/Data";
mkdir("$subDir") unless $fileCheck == 0;
my $subDir2 = cwd()."/Analysis";
mkdir("$subDir2") unless $fileCheck == 0;
my $genomeName = $ARGV[2];
my $trimMode = $ARGV[3];
my $trimLength = $ARGV[4];
my (%list,$chr,$min,$max);
if (defined $ARGV[8]) {
	open my $IN, '<', $ARGV[8] or die "$!";
	<$IN> for (1..1);
	while (<$IN>) {
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
	open my $IN2, '<', $file or die "$!";	#Open and read input .SAM file(s)
	(my $strain = $file) =~ s/_[^_]+$//;	#Strain-name
	(my $mode = $ARGV[0]) =~ s/\.SAM//;		#Alignment-mode
	print "$strain\n";
	my $outfile = "Data.".$genomeName."_".$strain."_".$mode.$outExt;	#Output files
	my $outfile2 = "Data.".$genomeName."_".$strain."_Ambiguous".$outExt;
	my ($OUT, $OUT2, $OUT3, $OUT4);
	open $OUT, '>', "$subDir/$outfile" or die "$!";
	open $OUT2, '>>', "$subDir/$outfile2" or die "$!";
	print $OUT "PairID\tStrand\tChr\tPos\tReadLength\tCIGAR\tAdjustment\tRead-Flag\tmSize\tRepeat\n";
	print $OUT2 "PairID\tStrand\tChr\tPos\tReadLength\tCIGAR\tAdjustment\tRead-Flag\tmSize\tMD-Tag\n";
	if ($alignRound eq "Global.SAM" && $trimMode eq "Y") {
		my $outfile3 = $strain.$ARGV[5]."_unmapped_trimmed.fastq";	#Unmapped R1 FASTQ file
		my $outfile4 = $strain.$ARGV[6]."_unmapped_trimmed.fastq";	#Unmapped R2 FASTQ file
		open $OUT3, '>', "$outfile3" or die "$!";
		open $OUT4, '>', "$outfile4" or die "$!";
	}
	my ($ID,$AID,$A,$B,$R1var,$R2var,$Wflag,$Cflag,%Size,%DoubleCut,$multimap);
	my $repeat = (0);
	while (<$IN2>) {	#For-each .SAM record
		chomp $_;
		next if /^\s*@/; #Skip .SAM headerlines
		my @F = split("\t", $_);	#Split each tab-delimited field
		my $orientation = $F[3]-$F[7];	#Discard atypical read-orientations
		if ($F[1] == 99 && $orientation > 0 || $F[1] == 83 && $orientation < 0) {
			my $skipline = <$IN2>;
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
			my $partner = <$IN2>;
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
					$DoubleCut{$chr}{$pos}{$posp}{"99"}++;
				}
			} elsif ($F[1] == 83) {
				$Wflag = 163; $Cflag = 83;
				if ($ARGV[0] eq "Global.SAM") {
					$Size{$mlength}++;
					$DoubleCut{$chr}{$pos}{$posp}{"163"}++;
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
				print $OUT join("\t",$ID,"w",$chr,$pos,$rl,$cigar,0-$Lclip,$Wflag,$mlength,$rflag),"\n",join("\t",$ID,"c",$chrp,$posp,$rlp,$cigarp,$Rclip+$cc-1,$Cflag,$mlength,$rflag),"\n";
			}
		}
	}
	close $IN2;
	close $OUT;
	close $OUT2;
	my $log = cwd()."/Logs/Log.".$strain.$outExt;
	my ($TRP,$MP,$MM,$TM,$AR,$MMR,$CPR,$ARR);
	open my $IN3, '+<', $log or die "$!";
	while (<$IN3>) {	#For-each .SAM record
		chomp $_;
		if (/(\d+) \(\d+.\d+\%\) were paired/) { $TRP = $1; }
		if (/(\d+) \(\d+.\d+\%\) aligned concordantly exactly/) { $MP = $1; }
		if (/(\d+) \(\d+.\d+\%\) aligned concordantly >/) { $MM = $1; }
	}
	print $IN3 "\n---------------------------------------\nCALL STATS\n---------------------------------------\n";
	print $IN3 "Total Pairs:\t",$TRP,"\n";
	printf($IN3 "%s\t%d\n","Valid Pairs:",$ID);
	printf($IN3 "%s\t%d\n","Repeat Pairs:",$repeat);
	printf($IN3 "%s\t%d\n","Ambig Pairs:",$AID);
	if ($ARGV[0] eq "Global.SAM") {
		my $chartfile = cwd()."/Logs/"."Plot.".$strain.".eps";
		my $chart = Chart::Gnuplot->new(output => $chartfile, title => "$strain Library Stats", ylabel => 'Count', grid => 'on', "style histogram" => "rowstack", legend => {position => 'outside top right', order => "vertical", align => 'left', width => -3, height => 0.8, sample => {position => "left"}}, imagesize => "0.6, 0.775", bg => {color => "#FFFFFF"});
		my @x = qw(Alignment | Mapping | Calling);
		my @y = ($MM+$MP,0,0,0,0); my @y2 = ($TRP-($MM+$MP),0,0,0,0);
		my @y3 = (0,0,$MM,0,0); my @y4 = (0,0,$MP,0,0);
		my @y5 = (0,0,0,0,$AID); my @y6 = (0,0,0,0,$repeat); my @y7 = (0,0,0,0,$ID-$repeat);
		$AR = (int(($MP/$TRP)*10000)/100); $MMR = (int(($MM/$TRP)*10000)/100);
		$ARR = (int(($AID/($MP+$MM))*10000)/100); my $DR = (int(($repeat/($MP+$MM))*10000)/100);
		$CPR = (int(($ID/($MP+$MM))*10000)/100)-$DR;
		my $uniqueRate = (int(($MP/($MP+$MM))*10000)/100);
		my $multiRate = (int(($MM/($MP+$MM))*10000)/100);
		print $OUT5 join("\t",$strain,$TRP,$MP+$MM,$MMR+$AR,$MP,$uniqueRate,$MM,$multiRate,$ID,$CPR,$AID,$ARR,$repeat,$DR),"\n";
		my $dataset = Chart::Gnuplot::DataSet->new(xdata => \@x, ydata => \@y, fill => {density => 0.5}, title => "Total Mapped Pairs (M)", style => "histograms", color => '#0066CC', border => { color => 'white'});
		my $dataset2 = Chart::Gnuplot::DataSet->new(xdata => \@x, ydata => \@y2, fill => {density => 0.5}, title => "Total Read Pairs", style => "histograms", color => '#D3D3D3', border => { color => 'white'});
		my $dataset3 = Chart::Gnuplot::DataSet->new(xdata => \@x, ydata => \@y3, fill => {density => 0.5}, title => "Mapped Pairs (Multi) ($multiRate% of M)", style => "histograms", color => '#FFC119', border => { color => 'white'});
		my $dataset4 = Chart::Gnuplot::DataSet->new(xdata => \@x, ydata => \@y4, fill => {density => 0.5}, title => "Mapped Pairs (Unique) ($uniqueRate% of M)", style => "histograms", color => '#CC0000', border => { color => 'white'});
		my $dataset5 = Chart::Gnuplot::DataSet->new(xdata => \@x, ydata => \@y5, fill => {density => 0.5}, title => "Ambiguous Pairs ($ARR% of M)", style => "histograms", color => '#8E44AD', border => { color => 'white'});
		my $dataset6 = Chart::Gnuplot::DataSet->new(xdata => \@x, ydata => \@y7, fill => {density => 0.5}, title => "Valid Pairs ($CPR% of M)", style => "histograms", color => '#17A817', border => { color => 'white'});
		if ($repeat >0) {
			my $dataset7 = Chart::Gnuplot::DataSet->new(xdata => \@x, ydata => \@y6, fill => {density => 0.4}, title => "Repeat Pairs ($DR% of M)", style => "histograms", color => '#000000', border => { color => 'white'});
			$chart->plot2d($dataset,$dataset2,$dataset3,$dataset4,$dataset5,$dataset7,$dataset6);
		} else {
			$chart->plot2d($dataset,$dataset2,$dataset3,$dataset4,$dataset5,$dataset6);
		}
		my $outfile6 = "Sizes.".$genomeName."_".$strain.$outExt;
		my $outfile7 = "Coordinates.".$genomeName."_".$strain.$outExt;
		open my $OUT6, '>', "$subDir2/$outfile6" or die "$!";
		open my $OUT7, '>', "$subDir2/$outfile7" or die "$!";
		print $OUT6 "Size\tFreq\n";
		print $OUT7 "Chr\tCoord-A\tCoord-B\tR1W Freq\tR1C Freq\n";
		foreach my $key (sort {$a <=> $b} keys %Size) {
			print $OUT6 "$key\t$Size{$key}\n";
		}
		foreach my $chrn (nsort keys %DoubleCut) {
			foreach my $coord1 (sort {$a <=> $b} keys %{$DoubleCut{$chrn}}) {
				foreach my $coord2 (sort {$a <=> $b} keys %{$DoubleCut{$chrn}{$coord1}}) {
						print $OUT7 join("\t",$chrn, $coord1,$coord2,map $_ // 0, @{$DoubleCut{$chrn}{$coord1}{$coord2}}{qw{99 163}});
						print $OUT7 "\n";
				}
			}
		}
	}
	if ($ARGV[1] eq "Y") {
		chomp($file);
		unlink($file);
	}
}
