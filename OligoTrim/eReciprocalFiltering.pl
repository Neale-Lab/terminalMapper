use strict;
use warnings;
use Cwd;
my $file = $ARGV[0];
my $outExt = '.txt';
open my $IN, '<', $file or die "$!";
my %expanded;
my $outfile = "Expanded_".$file;
open my $OUT, '>>', $outfile or die "$!";
my $header = <$IN>;
while (<$IN>) {
	chomp $_;
  my $splitR1W; my $splitR1C;
  my @F = split("\t", $_);
  if ($F[3]==0) {
    $splitR1W = 0
  } else {
    $splitR1W = $F[3]/9;
  }
  if ($F[4]==0) {
    $splitR1C = 0
  } else {
    $splitR1C = $F[4]/9;
  }
	$expanded{$F[0]}{$F[1]}{$F[2]}{"R1C"} += $splitR1C;
	$expanded{$F[0]}{$F[1]}{$F[2]+1}{"R1C"} += $splitR1C;
	$expanded{$F[0]}{$F[1]}{$F[2]-1}{"R1C"} += $splitR1C;
	$expanded{$F[0]}{$F[1]-1}{$F[2]}{"R1C"} += $splitR1C;
	$expanded{$F[0]}{$F[1]-1}{$F[2]+1}{"R1C"} += $splitR1C;
	$expanded{$F[0]}{$F[1]-1}{$F[2]-1}{"R1C"} += $splitR1C;
	$expanded{$F[0]}{$F[1]+1}{$F[2]}{"R1C"} += $splitR1C;
	$expanded{$F[0]}{$F[1]+1}{$F[2]+1}{"R1C"} += $splitR1C;
	$expanded{$F[0]}{$F[1]+1}{$F[2]-1}{"R1C"} += $splitR1C;
	$expanded{$F[0]}{$F[1]}{$F[2]}{"R1W"} += $splitR1W;
	$expanded{$F[0]}{$F[1]}{$F[2]+1}{"R1W"} += $splitR1W;
	$expanded{$F[0]}{$F[1]}{$F[2]-1}{"R1W"} += $splitR1W;
	$expanded{$F[0]}{$F[1]-1}{$F[2]}{"R1W"} += $splitR1W;
	$expanded{$F[0]}{$F[1]-1}{$F[2]+1}{"R1W"} += $splitR1W;
	$expanded{$F[0]}{$F[1]-1}{$F[2]-1}{"R1W"} += $splitR1W;
	$expanded{$F[0]}{$F[1]+1}{$F[2]}{"R1W"} += $splitR1W;
	$expanded{$F[0]}{$F[1]+1}{$F[2]+1}{"R1W"} += $splitR1W;
	$expanded{$F[0]}{$F[1]+1}{$F[2]-1}{"R1W"} += $splitR1W;
}
foreach my $chr (sort keys %expanded) {
	foreach my $coordA (sort {$a<=>$b} keys %{$expanded{$chr}}) {
  	foreach my $coordB (sort {$a<=>$b} keys %{$expanded{$chr}{$coordA}}) {
    	if (defined $expanded{$chr}{$coordA}{$coordB}{'R1W'} and defined $expanded{$chr}{$coordA}{$coordB}{'R1C'}) {
				if ($expanded{$chr}{$coordA}{$coordB}{'R1W'} > 0 and $expanded{$chr}{$coordA}{$coordB}{'R1C'} > 0) {
    		print $OUT "$chr\t$coordA\t$coordB\t$expanded{$chr}{$coordA}{$coordB}{'R1W'}\t$expanded{$chr}{$coordA}{$coordB}{'R1C'}\n";
			}
  		}
  	}
	}
}
