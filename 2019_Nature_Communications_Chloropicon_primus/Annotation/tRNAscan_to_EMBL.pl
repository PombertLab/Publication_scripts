#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl tRNAscan_to_EMBL.pl tRNAscanOutput.txt';

die $usage unless @ARGV;

my %aa_3_to_1_letter_code = (
'Ala' => 'A','Asx' => 'B','Cys' => 'C','Asp' => 'D',
'Glu' => 'E','Phe' => 'F','Gly' => 'G','His' => 'H',
'Ile' => 'I','Lys' => 'K','Leu' => 'L','Met' => 'M',
'Asn' => 'N','Pyl' => 'O','Pro' => 'P','Gln' => 'Q',
'Arg' => 'R','Ser' => 'S','Thr' => 'T','seC' => 'U',
'Val' => 'V','Trp' => 'W','Xaa' => 'X','Tyr' => 'Y',
'Glx' => 'Z'
);

while (my $tRNAs = shift @ARGV) {
	chomp $tRNAs;
	open IN, "<$tRNAs";
	
	my $chromosome = undef;

	while (my $line = <IN>) {
		chomp $line;
		if ($line =~ /^(\w+)\s+\d+\s+(\d+)\s+(\d+)\s+(\w{3,6})\t(\w{3})\t\d+\t\d+\t(\d+.\d+)/){
			$chromosome = $1;
			my $start = $2;
			my $end = $3;
			my $aa = $4;
			my $AC = $5;
			my $cove = $6;
			$AC =~ s/T/U/g;
			my $ac = lc($AC);
			
			open OUT, ">>$chromosome.tRNAscan.embl";
			
			if ($start < $end) {
				print OUT "FT   tRNA             $start..$end\n";
				print OUT "FT                   /gene=trn$aa_3_to_1_letter_code{$aa}"."($ac)"."\n";
				print OUT "FT                   /product=\"tRNA-$aa($AC)\"\n";
			}
			elsif ($start > $end) {
				print OUT "FT   tRNA             complement($start..$end)\n";
				print OUT "FT                   /gene=trn$aa_3_to_1_letter_code{$aa}"."($ac)"."\n";
				print OUT "FT                   /product=\"tRNA-$aa($AC)\"\n";
			}
		}
	}
}
close IN;
close OUT;
exit;