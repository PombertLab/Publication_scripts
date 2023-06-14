#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2023

use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(GetOptions);
use Statistics::Lite qw(mean stddev);

=begin comment
## All bash command lines are in GCeq.sh
=end comment
=cut 

my $name = 'gceq.pl';
my $version = '0.2';
my $updated = '2023-05-19';

my $usage = <<"HOWTO";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Calculates GCeq values from VCF files based on equations from:
		# Evidence That Mutation Is Universally Biased towards AT in Bacteria
		# Hershberg R, Petrov DA.2010. PLoS Genetics. 6(9): e1001115.
		# https://doi.org/10.1371/journal.pgen.1001115

COMMAND		${name} \\
		  -f *.fasta \\
		  -v *.vcf \\
		  -o outfile \\
		  -verbose

OPTIONS:
-f (--fasta)	FASTA files used for reference (one will be chosen at random GCeq calculations)
-v (--vcf)	VCF files
-o (--out)	Output file name
-verbose	Adds verbosity
HOWTO
die "\n$usage\n" unless @ARGV;

my @fasta;
my @vcf;
my $verbose;
my $outfile;
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'v|vcf=s@{1,}' => \@vcf,
	'o|out=s' => \$outfile,
	'verbose' => \$verbose
);

######### Choosing a FASTA file at random
my %sequences;
my @prefixes;

foreach my $fasta (@fasta){
	
	open FASTA, '<', $fasta or die $!;

	my $basename = fileparse($fasta);
	my ($prefix) = $basename =~ /^(\w+)/;
	push (@prefixes, $prefix);
	
	while (my $line = <FASTA>){

		chomp $line;
		
		unless ($line =~ /^>/){
		
			## Concatenating all contigs/chromosomes
			## Not a problem for GCeq
			$sequences{$prefix}{'sequence'} .= $line;
		
		}

	}

}


## Number of AT and GC bases in sequences
foreach my $fasta (keys %sequences){

	my $sequence = uc($sequences{$fasta}{'sequence'});
	my $seqlen = length($sequence);

	my $g = $sequence =~ tr/G//;
	my $c = $sequence =~ tr/C//;
	my $a = $sequence =~ tr/A//;
	my $t = $sequence =~ tr/T//;

	my $gc_sites = $g + $c;
	my $at_sites = $a + $t;
	my $gc_percent = sprintf("%.2f", ($gc_sites/$seqlen) * 100);

	$sequences{$fasta}{'gc_sites'} = $gc_sites;
	$sequences{$fasta}{'at_sites'} = $at_sites;
	$sequences{$fasta}{'gc_percent'} = $gc_percent;

	## Making sure that the math is correct
	my $sub = $seqlen - $gc_sites - $at_sites;
	if ($verbose){
		print "\nMetrics for $fasta:\n";
		print "Total length: $seqlen"."\n";
		print "# of GC bases: $gc_sites"."\n";
		print "# of AT bases: $at_sites"."\n";
		print "Non-ATGC bases: $sub"."\n";
		print "\% GC: $gc_percent\n";
	}

}

## Random selection
my $rand_fasta = $prefixes[int rand (@prefixes)];
print "\n".'FASTA sequence randomly selected = '.$rand_fasta."\n\n";
my $sequence = $sequences{$rand_fasta}{'sequence'};


######################## Working on VCF file(s)

## GCeq equations from:
## Evidence That Mutation Is Universally Biased towards AT in Bacteria
## Ruth Hershberg, Dmitri A. Petrov
## 2010. PLoS Genetics
## https://doi.org/10.1371/journal.pgen.1001115

open OUT, '>', $outfile or die $!;
open MET, '>', "$outfile.metrics" or die $!;
print OUT "# VCF\tGCeq\n";
print MET "# Organism\tGC%\tGCeq (mean)\tGCeq (SD)\n";

my %metrics;

while (my $vcf = shift @vcf){

	## Removing self hits (datasets against same genome)
	## Not useful for GCeq; skews results
	my $basename = fileparse($vcf);
	my ($prefix) = $basename =~ /^(\w+)/;
	
	my @s = split (/$prefix\./,$basename);
	my $count = scalar @s - 1;

	next if ($count > 1);
		
	## Number of SNPs in VCF files GC -> AT and AT -> GC
	my $gc_at_snps = 0;
	my $at_gc_snps = 0;

	open VCF, "<", $vcf or die "$vcf: $!";

	while (my $line = <VCF>){
		
		if ($line =~ /^#/){
			next;
		}

		else {

			my @data = split("\t", $line);
			my $ref = $data[3];
			my $alt = $data[4];

			## skip indels (if any)
			if ( (length($ref) > 1) or (length($alt) > 1) ){
				next;
			}

			else {
				if ( ($ref eq 'G') or ($ref eq 'C') ){
					if ( ($alt eq 'A') or ($alt eq 'T')){
						$gc_at_snps++;
					}
				}
				if ( ($ref eq 'A') or ($ref eq 'T') ){
					if ( ($alt eq 'G') or ($alt eq 'C')){
						$at_gc_snps++;
					}
				}
			}
		
		}

	}


	## Ratios
	my $r_gc_at = $gc_at_snps / $sequences{$rand_fasta}{'gc_sites'};
	my $r_at_gc = $at_gc_snps / $sequences{$rand_fasta}{'at_sites'};

	## GCeq
	my $gceq;
	if ( ($r_at_gc == 0) and ($r_gc_at == 0) ){
		$gceq = 'N/A';
	} 
	else {
		$gceq = sprintf("%.3f", ($r_at_gc / ($r_at_gc + $r_gc_at))*100);
	}

	print OUT $vcf."\t".$gceq." %\n";
	push (@{$metrics{$prefix}{'gceq'}}, $gceq); 

}

## Calculating GCeq means and standard deviations
foreach my $fasta (sort (keys %metrics)){

	my $gc_percent = $sequences{$fasta}{'gc_percent'};
	my $mean  = mean(@{$metrics{$fasta}{'gceq'}});
	my $stdev  = stddev(@{$metrics{$fasta}{'gceq'}});

	print MET $fasta."\t".$gc_percent."\t".$mean."\t".$stdev."\n";

}

