#!/usr/bin/perl
# Pombert lab 2019
my $version = '0.2';
my $name = 'split_VCF.pl';
my $updated = '2023-04-24';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $options = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Splits VCF file(s) per contig/chromosome, excluding positions masked (N) by RepeatMasker in the reference.

COMMAND		${name} -v *.vcf -o ./

OPTIONS:
-v (--vcf)	VCF file(s) to parse
-o (--out)	Output directory [Default: ./]
OPTIONS
die "\n$options\n" unless @ARGV;

my @VCFs;
my $outdir = './';
GetOptions(
	'v|vcf=s@{1,}' => \@VCFs,
	'o|out=s' => \$outdir
);

unless (-d $outdir){
	mkdir ($outdir, 0755) or die $!;
}

## Parsing VCF files
while (my $vcf = shift@VCFs){

	open VCF, '<', $vcf or die $!;

	my $previous_contig = undef;
	my $filename;

	while (my $line = <VCF>){

		if ($line =~ /^#/){
			next;
		}
	
		chomp $line;
		my @columns = split("\t", $line);
		my $contig = $columns[0];
		my $ref_nt = $columns[3];
	
		if (!defined $previous_contig){

			$previous_contig = $contig;
			$filename = $outdir.'/'.$contig.'.vcf';

			## Using an autovivified anonymous file handle (referenced in $fh)
			open my $fh, '>', $filename or die $!;

			## Checking if position has been masked by RepeatMasker, if not continue...
			unless ($ref_nt eq 'N'){
				print $fh "$line\n";
			}

		}

		elsif ($contig ne $previous_contig) {

			$previous_contig = $contig;
			$filename = $outdir.'/'.$contig.'.vcf';

			open my $fh, '>', $filename or die $!;

			unless ($ref_nt eq 'N'){
				print $fh "$line\n";
			}

		}

		elsif ($contig eq $previous_contig){

			$previous_contig = $contig;
			$filename = $outdir.'/'.$contig.'.vcf';

			open my $fh, '>>', $filename or die $!;

			unless ($ref_nt eq 'N'){
				print $fh "$line\n";
			}

		}

	}

}