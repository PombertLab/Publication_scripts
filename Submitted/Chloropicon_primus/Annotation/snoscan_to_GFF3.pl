#!/usr/bin/perl

use strict;
use warnings;

my %snorna = ();
my $snow = 0;
	
while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.gff3";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>>\s+(\S+)\s+(\S+)\s+\((\d+)-(\d+)\).*\(([-+])/){
			my $target = $1;
			my $score = $2;
			my $start = $3;
			my $end = $4; 
			my $strand = $5;
			my $key = $target; $key .= '-'; $key .= $start; $key .= '-'; $key .= $end;

			if (exists $snorna{$key}){next;}
			else{
			$snow++;
			print OUT "$target"."\t"."SNOSCAN"."\t"."gene"."\t"."$start"."\t"."$end"."\t"."$score"."\t"."$strand"."\t".'.'."\t"."ID=snoRNA$snow".';'."Name=snoRNA$snow"."\n";
			print OUT "$target"."\t"."SNOSCAN"."\t"."snoRNA"."\t"."$start"."\t"."$end"."\t"."$score"."\t"."$strand"."\t".'.'."\t"."gene_id=snoRNA$snow".';'."Parent=snoRNA$snow".';'."transcript_id=$snow.t1"."\n";
			}
			$snorna{$key}= $key;
		}
	}
}