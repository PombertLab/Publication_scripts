#!/usr/bin/perl
# Pombert lab 2019
my $version = '0.2b';
my $name = 'sort_SNPs.pl';
my $updated = '2023-05-05';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;

## Defining options
my $options = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Parse allelic distributions in VCF files for easy plotting with R or Excel

COMMAND		sort_SNPs.pl \\
		  -min 10 \\
		  -max 90 \\
		  -vcf *.vcf \\
		  -noindel \\
		  -out ./

OPTIONS:
-min	Mininum variant allelic frequency to keep [Default: 10]
-max	Maximum variant allelic frequency to keep [Default: 90]
-vcf	VCF file(s) to parse
-out	Output directory [Default: ./]
-noindel	Skip indels		## Nanopore data generates too many indels
-extra	Creates additonal (extra) outputfiles (unsorted and gaussian)
OPTIONS
die "$options\n" unless @ARGV;

my $min = 10;
my $max = 90;
my @VCF;
my $outdir = './';
my $noindel;
my $extra;
GetOptions(
	'min=i' => \$min,
	'max=i' => \$max,
	'v|vcf=s@{1,}' => \@VCF,
	'out=s' => \$outdir,
	'noindel' => \$noindel,
	'extra' => \$extra
);

## Creating outdir if needed
unless (-d $outdir){
	mkdir ($outdir, 0755) or die $!;
}

## Parsing VCF file(s)
while (my $file = shift @VCF){

	open IN, '<', $file or die $!;
	my ($basename, $path) = fileparse($file);
	$basename =~ s/.vcf$//;

	open VCF, '>', "$outdir/$basename.sorted.$min.$max.vcf" or die $!;		## VCF parsed
	open RSO, '>', "$outdir/$basename.sorted.$min.$max.tsv" or die $!;		## File to graph with excel or R
	
	if ($extra){
		open RUN, '>', "$outdir/$basename.unsorted.$min.$max.tsv" or die $!; 	## File to graph with excel or R
		open GAUS, '>', "$outdir/$basename.gaussian.$min.$max.tsv" or die $!; 	## Data for Gaussian plot 
	} 

	my @percents;
	my %percents;
	
	while (my $line = <IN>){

		chomp $line;

		if ($line =~ /^#/){
			next;
		}

		else {

			my @columns = split("\t", $line);
			my $locus = $columns[0];
			my $position = $columns[1];
			my $id = $columns[2];
			my $ref = $columns[3];
			my $alt = $columns[4];
			my $format = $columns[9];

			my @stats = split(":", $format);
			my $allele_freq = $stats[6];
			
			if ($allele_freq){
				$allele_freq =~ s/%//;
				$allele_freq = sprintf("%.03d", $allele_freq);
			}
			
			## Skipping masked nucleotides
			if ($ref eq 'N'){
				next;
			}

			## Skipping indels
			if ($noindel){

				## Short deletions are really common in nanopore reads. 
				## Creates too many fake indels; homopolymers?
				if ($alt eq '.'){
					next;
				}

				if (length ($ref) != length ($alt)){
					next;
				} 
			}
			
			if ( ($allele_freq >= $min) && ($allele_freq <= $max) ){

				## Calculating frequency for alternate nucleotide, e.g.:
				#
				# diploid: 1 - 0.5 = 0.5
				# triploid: 1 - 0.66 = 0.33 (or vice versa)
				# tetraploid: 1 - 0.25 = 0.75 (or vice versa)
				#
				# Note, this only works if there is no contamination in the mix
				my $alternate = 100 - $allele_freq;

				push (@percents, $allele_freq);
				push (@percents, $alternate);
				$percents{$allele_freq} += 1;
				$percents{$alternate} += 1;

				print VCF "$line\n";

			} 
		}
	}
	
	my @unsort = @percents;
	@percents = sort @percents;
	@percents = reverse @percents;
	
	while (my $pc = shift @percents){
		print RSO "$pc\n";
	}
	
	if ($extra){
		while (my $un = shift @unsort){
			print RUN "$un\n";
		}

		for (keys %percents){
			print GAUS "$_\t$percents{$_}\n";
		}
	} 

}