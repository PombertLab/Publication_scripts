#!/usr/bin/perl

## Loads a list of ECUs from a text file, one ECU per line
## Then each corresponding ECI, ECII and ECIII sequence is loaded and put into their own arrays with one codon per line.
## Each aligned codon is compared between the three genes, and written to the appropriate output files.
## This is possible because the genes have been aligned. A single frame shift will kill this script.
## Note 1: that if alignments are not based on codons, each gap containing gene will need to be verified before running this script.
## Note 2: The sequence must be on one line without a fasta header. Use fasta_to_string.pl script to convert to strings.


use strict;
use warnings;

my $usage = 'perl Codon_SNPs.pl list.txt';

die $usage unless @ARGV;

while (my $txt = shift @ARGV) {

	open IN, "<$txt" or die "cannot open $txt";
	
		while (my $line = <IN>) {

			chomp $line;
		
			if ($line =~ /^(ECU\d+.*)/) {
			
				my $ECU = $1;
		
				my $seq1= "ECI_$ECU.string";	## match the file extention to that of the input file used
				my $seq2 = "ECII_$ECU.string";  ## match the file extention to that of the input file used
				my $seq3 = "ECIII_$ECU.string"; ## match the file extention to that of the input file used

				open (IN1, $seq1);
				open (IN2, $seq2);
				open (IN3, $seq3);

				open INVAR, ">$ECU.codon.invar";
				open SNPs, ">$ECU.codon.snps";
				open GAPs, ">$ECU.codon.gaps";

				my $strain1 = <IN1>;
				my $strain2 = <IN2>;
				my $strain3 = <IN3>;

				my @ECI = unpack ("(A3)*", $strain1);
				my @ECII = unpack ("(A3)*", $strain2);
				my @ECIII = unpack ("(A3)*", $strain3);

				my $length = scalar @ECI;
				my $end = $length-1;

				my $codon = 1;

				foreach (0..$end) {
	
					if (($ECI[$_] =~ /-/) || ($ECII[$_] =~ /-/) || ($ECIII[$_] =~ /-/)) {
					print GAPs "$codon\t$length\t$ECI[$_]\t$ECII[$_]\t$ECIII[$_]\n";
					$codon++; 
					}
				
					elsif (($ECI[$_] eq $ECII[$_]) && ($ECI[$_] eq $ECIII[$_])) {
					print INVAR "$codon\t$length\t$ECI[$_]\t$ECII[$_]\t$ECIII[$_]\n";
					$codon++; 
					}
		
					else {
					print SNPs "$codon\t$length\t$ECI[$_]\t$ECII[$_]\t$ECIII[$_]\n";
					$codon++; 
		
				}
			}	
		}

	close IN1;
	close IN2;
	close IN3;
	close INVAR;
	close SNPs;
	close GAPs;
	
	}

}

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)