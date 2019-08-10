#!/usr/bin/perl
## Pombert Lab, 2019
my $name = 'extract_PBSV_variants.pl';
my $version = 0.1;

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Extract stuctural variants from PBSV VCF files
EXAMPLE		extract_PBSV_variants.pl -m 500 -v *.vcf 

OPTIONS:
-m (--min)	Minimum size of the variant [default: 500]
-v (--vcf)	VCF file(s) to parse

OPTIONS
die $usage unless @ARGV;

my $min = 500;
my @vcf;
GetOptions(
	'm|min=i' => \$min,
	'v|vcf=s@{1,}' => \@vcf
);

while (my $vcf = shift@vcf){
	open VCF, "<$vcf" or die "\nCan't open VCF file named $vcf. Please check your command line...\n\n";
	$vcf =~ s/.vcf$//;
	open FASTA, ">$vcf.m$min.fasta";
	while (my $line = <VCF>){
		chomp $line;
		if ($line =~ /^#/){next;} ## Skipping headers
		else{
			my @columns = split("\t", $line); #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
			my $contig = $columns[0];
			my $variant = $columns[2];
			my $ref = $columns[3];
			my $alt = $columns[4];
			my $filter = $columns[6];
			if ($variant =~ /DEL/){
				my $len = length($ref);
				if ($len >= $min){
					print "Deletion $variant on $contig is $len bp\n";
					print FASTA ">$variant [$contig; DEL]\n";
					my @seq = unpack ("(A60)*", $ref);
					while (my $tmp = shift@seq){print FASTA "$tmp\n";}
				}
			}
			elsif ($variant =~ /INS/){
				my $len = length($alt);
				if ($len >= $min){
					print "Insertion $variant on $contig is $len bp\n";
					print FASTA ">$variant [$contig; INS]\n";
					my @seq = unpack ("(A60)*", $alt);
					while (my $tmp = shift@seq){print FASTA "$tmp\n";}
				}
			}
		}
	}
	close VCF;
	close FASTA;
}