#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl Aragorn_to_EMBL.pl AragornOutput.txt';

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
		if ($line =~ />(Chromosome_\d+)/){
		$chromosome = $1;
		}
		if ($line =~ /^\d+\s+tRNA-(\w{3}).*c\[(\d+),(\d+)\]\s+\d+.*\((\w{3})\)/){
			my $aa = $1;
			my $start = $2;
			my $end = $3;
			my $ac = $4;
			$ac =~ s/t/u/g;
			my $AC = uc($ac);

			open OUT, ">>$chromosome.aragorn.embl";
			print OUT "FT   tRNA             complement($start..$end)\n";
			print OUT "FT                   /gene=trn$aa_3_to_1_letter_code{$aa}"."($ac)"."\n";
			print OUT "FT                   /product=\"tRNA-$aa($AC)\"\n";
			print OUT "FT                   /note=\"inferred with Aragorn\"\n";
		}
		elsif ($line =~ /^\d+\s+tRNA-(\w{3}).*\[(\d+),(\d+)\]\s+\d+.*\((\w{3})\)/){
			my $aa = $1;
			my $start = $2;
			my $end = $3;
			my $ac = $4;
			$ac =~ s/t/u/g;
			my $AC = uc($ac);			

			open OUT, ">>$chromosome.aragorn.embl";
			print OUT "FT   tRNA             $start..$end\n";
			print OUT "FT                   /gene=trn$aa_3_to_1_letter_code{$aa}"."($ac)"."\n";
			print OUT "FT                   /product=\"tRNA-$aa($AC)\"\n";
			print OUT "FT                   /note=\"inferred with Aragorn\"\n";
		}
	}
}

close IN;
close OUT;
exit;