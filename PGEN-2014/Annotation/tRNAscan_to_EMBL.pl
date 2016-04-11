#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl tRNAscan_to_EMBL.pl tRNAscanOutput.txt';

die $usage unless @ARGV;

my %aa_3_to_1_letter_code = (
'Ala' => 'A','Asx' => 'B','Cys' => 'C','Asp' => 'D',
'Glu' => 'E','Phe' => 'F','Gly' => 'G','His' => 'H',
'Ile' => 'I','Lys' => 'K','Leu' => 'L','Met' => 'M',
'Asn' => 'N','Pro' => 'P','Gln' => 'Q','Arg' => 'R',
'Ser' => 'S','Thr' => 'T','Val' => 'V','Trp' => 'W',
'Xaa' => 'X','Tyr' => 'Y','Glx' => 'Z'
);

my %anticodon = (
'CGA' => 'cga','CGG' => 'cgg','CGT' => 'cgu','CGG' => 'cgg',
'GCA' => 'gca','GCG' => 'gcg','GCT' => 'gcu','GCC' => 'gcc',
'TCT' => 'ucu','TCC' => 'ucc','TTA' => 'uua','TTG' => 'uug',
'CTA' => 'cua','CTG' => 'cug','ACA' => 'aca','ACG' => 'acg',
'CTT' => 'cuu','CTC' => 'cuc','GTT' => 'guu','GTC' => 'guc',
'CCA' => 'cca','CCG' => 'ccg','CCC' => 'ccc','CCT' => 'ccu',
'GTA' => 'gua','GTG' => 'gug','TAA' => 'uaa','TAG' => 'uag',
'TAT' => 'uau','AAT' => 'aau','AAC' => 'aac','GAA' => 'gaa',
'GAG' => 'gag','GAT' => 'gau','GAC' => 'gac','TTT' => 'uuu',
'TTC' => 'uuc','TAC' => 'uac','AAA' => 'aaa','AAG' => 'aag',
'GGA' => 'gga','GGG' => 'ggg','GGT' => 'ggu','GGC' => 'ggc',
'AGA' => 'aga','AGG' => 'agg','AGT' => 'agu','AGC' => 'agc',
'TCA' => 'uca','TCG' => 'ucg','ATG' => 'end','ATT' => 'end',
'ACT' => 'end','TGA' => 'uga','TGG' => 'ugg','TGT' => 'ugu',
'TGC' => 'ugc','ACC' => 'acc','ATA' => 'aua','ATG' => 'aug',
'CAA' => 'caa','CAG' => 'cag','CAT' => 'cau','CAC' => 'cac'
);

while (my $tRNAs = shift @ARGV) {
	chomp $tRNAs;
	open IN, "<$tRNAs";
	open OUT, ">$tRNAs.embl";

	while (my $line = <IN>) {
		chomp $line;
		if ($line =~ /^\w+\s+\d+\s+(\d+)\s+(\d+)\s+(\w{3})\t(\w{3})\t\d\t\d\t(\d+.\d+)/){
			my $start = $1;
			my $end = $2;
			my $aa = $3;
			my $ac = $4;
			my $cove = $5;
			
			if ($start < $end) {
				print OUT "FT   tRNA             $start..$end\n";
				print OUT "FT                   /gene=trn$aa_3_to_1_letter_code{$aa}"."($anticodon{$ac})"."\n";
				print OUT "FT                   /product=\"tRNA-$aa($ac)\"\n";
			}
			elsif ($start > $end) {
				print OUT "FT   tRNA             complement($start..$end)\n";
				print OUT "FT                   /gene=trn$aa_3_to_1_letter_code{$aa}"."($anticodon{$ac})"."\n";
				print OUT "FT                   /product=\"tRNA-$aa($ac)\"\n";
			}
		}
	}
}
close IN;
close OUT;
exit;