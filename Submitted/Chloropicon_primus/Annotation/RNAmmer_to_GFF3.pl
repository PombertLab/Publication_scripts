#!/usr/bin/perl
## Converts the output of RNAmmer GFF2 output to the proper GFF3 format for loading into webApollo

use strict;
use warnings;

my $usage = 'USAGE = perl RNAmmer_to_GFF3.pl *.gff2';
die "$usage\n" unless @ARGV;

my $rRNA = 0;

while (my $file = shift@ARGV){
	open RRNA, "<$file";
	$file =~ s/.gff2$//;
	open OUT, ">$file.gff3";
		
	## Converting rRNAs to GFF3
	while (my $line = <RRNA>){
		chomp $line;
		if ($line =~ /^#/){next;} ## disregard comments
		elsif ($line =~ /^(\S+)\s+\S+\s+\S+\s(\d+)\s+(\d+)\s+(\S+)\s+([+-])\s+\.\s+(\S+)/){
			my $target = $1;
			my $start = $2;
			my $end = $3;
			my $score = $4;
			my $strand = $5;
			my $gene = $6;
			$rRNA++;

			if ($strand eq '+'){
				print OUT "$target"."\t"."RNAMMER"."\t"."rRNA"."\t"."$start"."\t"."$end"."\t"."$score"."\t".'+'."\t".'.'."\t"."ID=rRNA$rRNA".';'."Name=rRNA$rRNA".';'."Note=$gene"."\n";
			}
			elsif ($strand eq '-'){
				print OUT "$target"."\t"."RNAMMER"."\t"."rRNA"."\t"."$start"."\t"."$end"."\t"."$score"."\t".'-'."\t".'.'."\t"."ID=rRNA$rRNA".';'."Name=rRNA$rRNA".';'."Note=$gene"."\n";
			}
		}
	}
}

