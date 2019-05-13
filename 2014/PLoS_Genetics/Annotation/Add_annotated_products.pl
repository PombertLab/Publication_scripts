#!/usr/bin/perl

use strict;
use warnings;

## Open *.tbl.3 files from my perl gauntlet, and then look up locus_tags in a set of defined text files (PFAMs.txt, Evalues.txt Products.txt)
## These files will be used as hashes

my $usage = 'perl Add_annotated_products.pl *.tbl3';

die $usage unless @ARGV;

open (PRODUCT, "<Products.txt") or die "Can't open Products!"; ## The hash input is a gene & product tab-delimited list, one per line,.
my %product = ();

open (EVALUE, "<Evalues.txt") or die "Can't open Evalues!"; ## The hash input is a gene & product tab-delimited list, one per line,.
my %evalue = ();

open (PFAM, "<PFAMs.txt") or die "Can't open PFAMs!"; ## The hash input is a gene & product tab-delimited list, one per line,.
my %pfam = ();

#### Fill the hashes with the genes (keys) and products (values).

	while (my $line = <PRODUCT>) {
		chomp $line;
		
		if ($line =~ /^(H632_c\d{1,}\.\d{1,})\t(.*$)/) {
			my $key = $1;
			my $value = $2;
			$product{"$key"} = "$value";
		}
	}
	
	while (my $line = <EVALUE>) {
		chomp $line;
		
		if ($line =~ /^(H632_c\d{1,}\.\d{1,})\t(.*$)/) {
			my $key = $1;
			my $value = $2;
			$evalue{"$key"} = "$value";
		}
	}
	
	while (my $line = <PFAM>) {
		chomp $line;
		
		if ($line =~ /^(H632_c\d{1,}\.\d{1,})\t(\S+)/) {
			my $key = $1;
			my $value = $2;
			$pfam{"$key"} = "$value";
		}
	}

#### Work on the query file. ###########

while (my $tbl = shift @ARGV) {

	open IN, "<$tbl" or die "cannot open $tbl";
	$tbl =~ s/\.tbl3$//;
	open OUT, ">$tbl.tbl4";

	while (my $line = <IN>)  {
		
		chomp $line;
		
		if ($line ~~ /^\t\t\tlocus_tag\t(H632_c\d{1,}\.\d{1,})/) {
			my $locus_tag = $1;
					
			if ($pfam{$locus_tag} =~ /PF/){
				print OUT "$line\n";
				print OUT "\t\t\tproduct\t$product{$locus_tag}\n";
				print OUT "\t\t\tnote\tsimilar to PFAM: $pfam{$locus_tag}, $evalue{$locus_tag}\n";
				#~ print OUT "\t\t\tnote\tPFAM: $pfam{$locus_tag}, E-value: $evalue{$locus_tag}\n";
			}else{
				print OUT "$line\n";
				print OUT "\t\t\tproduct\thypothetical protein\n";
			}
		}
		else {
			print OUT "$line\n";
		}
	}
		
	close IN;
}

exit;
