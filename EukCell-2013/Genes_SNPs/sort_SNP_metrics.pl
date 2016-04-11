#!/usr/bin/perl

## Bin SNPs from the cat_gene_SNPs.pl output according to desired thresholds.
## Note: To bin the cat_gene_SNPs.pl output by percentiles, simply open the file in a LibreOffice or Excel spreadsheet and sort accordingly. 

use strict;
use warnings;

my $usage = 'perl sort_SNP_metrics.pl all_gene_SNPs_metrics.txt';

die $usage unless @ARGV;

open ZERO, ">all_gene_SNPs_metrics.zero";
open LOW, ">all_gene_SNPs_metrics.low";
open MID, ">all_gene_SNPs_metrics.mid";
open HIGH, ">all_gene_SNPs_metrics.high";

my $low = 0;	## set here the desired ratios to sort the genes according to their relative amount of SNPs
my $mid = 0.001;	## set here the desired ratios to sort the genes according to their relative amount of SNPs
my $high = 0.01;	## set here the desired ratios to sort the genes according to their relative amount of SNPs

while (my $metrics = shift @ARGV) {

	open IN, "<$metrics" or die "cannot open $metrics";
	
		while (my $line = <IN>) {

			chomp $line;
		
			if ($line =~ /^(ECU\d+\S+)\t(\d+)\t(\d+)\t(\d+.*)/) {
				
				my $gene = $1;
				my $length = $2;
				my $SNPs = $3;
				my $ratio = $4;
			
				if ($SNPs == 0){
				print ZERO "$line\n";
				}
				elsif ($ratio >= $high){
				print HIGH "$line\n";
				}
				elsif (($ratio >= $mid) && ($ratio < $high)){
				print MID "$line\n";
				}
				elsif (($ratio > 0) && ($ratio < $mid)){
				print LOW "$line\n";
				}
			}
		}
	close IN;
}

close ZERO;
close LOW;
close MID;
close HIGH;

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)