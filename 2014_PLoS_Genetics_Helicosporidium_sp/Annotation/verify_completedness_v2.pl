#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl verify_completedness.pl *.feat';

die $usage unless @ARGV;

## verify the completion in 5' & 3' end of the genes predicted by the MAKER2 gauntlet based on the presence of start and stop codons ##
## starts from the .feat files generated with extractGFF_features.pl ##
## NOTE: the codons on the reverse strand are searched for by their reverse complement, much easier to script since no need to play with the strings ##

while (my $feat = shift @ARGV) {

	open IN, "<$feat" or die "cannot open $feat";
	$feat =~ s/\.feat$//;
	open OUT, ">$feat.ends";
	
	my @genes = ();

	while (my $line = <IN>) {

		chomp $line;

		if ($line =~ /^[+-]\tmatch\t(\d+)\t(\d+)\t(\S+)$/) {
			push (@genes, $line);
		}
		else {
			next;
		}
	}	
	close IN;
	
	open STRING, "<$feat.string" or die "cannot open $feat.string";
	my $string = <STRING>;
	
	foreach my $genes (@genes) {
		if ($genes =~ /^(\+)\tmatch\t(\d+)\t(\d+)\t(\S+)$/) {
			my $strand = $1;
			my $start = ($2-1);
			my $end = ($3-1);
			my $locus = $4;
			my $start_codon = substr($string, $start, 3);
			my $stop_codon = substr($string, ($end-2), 3);
			
				if ($start_codon eq 'ATG'&& ($stop_codon eq 'TAA' || $stop_codon eq 'TAG' || $stop_codon eq 'TGA')) { ## bug in precedence for stop_codon?
					print OUT "$strand\t$locus\tcomplete\t$start_codon\t$stop_codon\t$2\t$3\n";
				}
				elsif ($start_codon ne 'ATG' && ($stop_codon eq 'TAA' || $stop_codon eq 'TAG' || $stop_codon eq 'TGA')) {
					print OUT "$strand\t$locus\t5end_missing\t$start_codon\t$stop_codon\t$2\t$3\n";
				}
				elsif ($start_codon eq 'ATG' && ($stop_codon ne 'TAA' && $stop_codon ne 'TAG' && $stop_codon ne 'TGA')) {
					print OUT "$strand\t$locus\t3end_missing\t$start_codon\t$stop_codon\t$2\t$3\n";
				}
				elsif ($start_codon ne 'ATG' && ($stop_codon ne 'TAA' && $stop_codon ne 'TAG' && $stop_codon ne 'TGA')) {
					print OUT "$strand\t$locus\tboth_ends_missing\t$start_codon\t$stop_codon\t$2\t$3\n";
				}
		}
		elsif ($genes =~ /^(\-)\tmatch\t(\d+)\t(\d+)\t(\S+)$/) {
			my $strand = $1;
			my $start = ($3-1);
			my $end = ($2-1);
			my $locus = $4;
			my $start_codon = substr($string, ($start-2), 3);
			my $stop_codon = substr($string, $end, 3);
			
						
				if ($start_codon eq 'CAT' &&  ($stop_codon eq 'TCA' || $stop_codon eq 'CTA' || $stop_codon eq 'TTA')) { ## bug in precedence for stop_codon?
					print OUT "$strand\t$locus\tcomplete\t$start_codon\t$stop_codon\t$3\t$2\n";
				}
				elsif ($start_codon ne 'CAT' && ($stop_codon eq 'TCA' || $stop_codon eq 'CTA' || $stop_codon eq 'TTA')) {
					print OUT "$strand\t$locus\t5end_missing\t$start_codon\t$stop_codon\t$3\t$2\n";
				}
				elsif ($start_codon eq 'CAT' && ($stop_codon ne 'TCA' && $stop_codon ne 'CTA' && $stop_codon ne 'TTA')) {
					print OUT "$strand\t$locus\t3end_missing\t$start_codon\t$stop_codon\t$3\t$2\n";
				}
				elsif ($start_codon ne 'CAT' && ($stop_codon ne 'TCA' && $stop_codon ne 'CTA' && $stop_codon ne 'TTA')) {
					print OUT "$strand\t$locus\tboth_ends_missing\t$start_codon\t$stop_codon\t$3\t$2\n";
				}
		}
		
	}
	
	close STRING;
	
}

exit;
