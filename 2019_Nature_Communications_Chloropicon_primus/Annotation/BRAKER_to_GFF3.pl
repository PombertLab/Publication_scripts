#!/usr/bin/perl
## Converts the GFF output of the BRAKER1 pipeline to the proper GFF3 format for loading into webApollo

use strict;
use warnings;

my $usage = 'USAGE = perl BRAKER_to_GFF3.pl *.gff';
die "$usage\n" unless @ARGV;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.gff$//;
	open OUT, ">$file.formatted.gff";
	
	my $gene = undef;
	my $ts_count = 0;
	my $cds_count = 0;
	
	while (my $line = <IN>){
		chomp $line;
		## Discarding comments
		if ($line =~ /^\#/){next;}
		## Printing genes
		elsif ($line =~ /^(\S+)\t(AUGUSTUS)\t(gene)\t(\d+)\t(\d+)\t(\S+)\t([-+])\t(\S)\t(\S+)/){
			my $seqid = $1;	## http://www.sequenceontology.org/gff3.shtml
			my $source = $2;
			my $type = $3;
			my $start = $4;
			my $end = $5;
			my $score = $6;
			my $strand = $7;
			my $phase = $8;
			$gene = $9;
			print OUT "$seqid"."\t"."$source"."\t"."match"."\t"."$start"."\t"."$end"."\t"."$score"."\t"."$strand"."\t"."$phase"."\t"."ID=$gene".';'."Name=$gene"."\n";
			$ts_count = 0;
			$cds_count = 0;
		}
		# Printing transcripts
		# elsif ($line =~ /^(\S+)\t(AUGUSTUS)\t(transcript)\t(\d+)\t(\d+)\t(\S+)\t([-+])\t(\S)/){
			# $ts_count++;
			# my $seqid = $1;	## http://www.sequenceontology.org/gff3.shtml
			# my $source = $2;
			# my $type = $3;
			# my $start = $4;
			# my $end = $5;
			# my $score = $6;
			# my $strand = $7;
			# my $phase = $8;
			# print OUT "$seqid"."\t"."$source"."\t"."$type"."\t"."$start"."\t"."$end"."\t"."$score"."\t"."$strand"."\t"."$phase"."\t"."ID=transcript$ts_count".';'."Parent=$gene".';'."transcript_id=mRNA$gene"."\n";
		# }
		## Printing CDS
		elsif ($line =~ /^(\S+)\t(AUGUSTUS)\t(CDS)\t(\d+)\t(\d+)\t(\S+)\t([-+])\t(\S)/){
			$cds_count++;
			my $seqid = $1;	## http://www.sequenceontology.org/gff3.shtml
			my $source = $2;
			my $type = $3;
			my $start = $4;
			my $end = $5;
			my $score = $6;
			my $strand = $7;
			my $phase = $8;
			print OUT "$seqid"."\t"."$source"."\t"."match_part"."\t"."$start"."\t"."$end"."\t"."$score"."\t"."$strand"."\t"."$phase"."\t"."gene_id=$gene".';'."Parent=$gene".';'."transcript_id=$gene.t1"."\n";
		}
	}
}

