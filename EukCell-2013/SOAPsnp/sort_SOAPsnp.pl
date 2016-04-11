#!/usr/bin/perl

## This script sorts the output of SOAPsnp.

use strict;
use warnings;

my $usage = 'perl sort_SOAPsnp.pl *.soapsnp.out';

die $usage unless @ARGV;

######Enter the desired minimum quality score below #########

my $minQS = 30;	## 

######Enter the desired minimum quality score above #########

while (my $soapsnp = shift @ARGV) {

	open IN, "<$soapsnp" or die "cannot open $soapsnp file";
	$soapsnp =~ s/\.out$//;
	open VOID, ">$soapsnp.void";			## Positions that are not covered by any sequences are put in the void file (useful to detect lack of coverage and indels).
	open SNPQC, ">$soapsnp.snps.qc";		## SNPs that pass the quality filter
	open SNPNOQC, ">$soapsnp.snps.noqc";	## SNPs that don't pass the quality filter
	open INVAR, ">$soapsnp.invar";			## Invariables sites
		
	while (my $line = <IN>) {

		chomp $line;
		
		if ($line =~ /^(gi\S+)\t(\d+)\t(\w)\t(\w)\t(\d+)\t(\w)\t(\d+)\t(\d+)\t(\d+)\t(\w)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+.\d+)\t(\d+.\d+)\t(\d+)/) {
			
			my $chromosome = $1;		#Chromosome ID
			my $locus = $2;			#Coordinate on chromosome, start from 1
			my $ref_genome = $3;		#Reference genotype
			my $query = $4; 			#Consensus genotype
			my $queryQS	= $5;		#Quality score of consensus genotype
			my $bb = $6;				#Best base
			my $bbQS = $7;			#Average quality score of best base
			my $bbCount = $8; 		#Count of uniquely mapped best base
			my $bbAll = $9;			#Count of all mapped best base
			my $Sndbb = $10;			#Second best bases
			my $SndbbQS = $11;		#Average quality score of second best base
			my $SndbbCount = $12;		#Count of uniquely mapped second best base
			my $SndbbAll = $13;		#Count of all mapped second best base
			my $SeqDepth = $14;		#Sequencing depth of the site
			my $pvaluew = $15;		#Rank sum test p_value
			my $AvCpNum = $16;		#Average copy number of nearby region
			my $dbSNP = $17;			#Whether the site is a dbSNP
			
			if ($SeqDepth == 0) {
				print VOID "$line\n";
			}
			
			elsif (($SeqDepth > 0) && ($query ne $ref_genome) && ($bbQS >= $minQS)){
				print SNPQC "$line\n";
			}
			
			elsif (($SeqDepth > 0) && ($query ne $ref_genome) && ($bbQS < $minQS)){
				print SNPNOQC "$line\n";
			}

			elsif (($SeqDepth > 0) && ($query eq $ref_genome)){
				print INVAR "$line\n";
			}

		}
		
		else {
			next;
		}
	}
	close IN;
	close VOID;
	close SNPQC;
	close SNPNOQC;
	close INVAR;
}

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)