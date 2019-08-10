#!/usr/bin/perl
## Pombert Lab, 2017
my $name = 'VCF_to_Circos.pl';
my $version = 0.1;

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Generates variants distribution from VarScan2 VCF files using sliding windows to plot with Circos
EXAMPLE		VCF_to_Circos.pl -v *.vcf -f *.fasta

OPTIONS:
-f (--fasta)	Input files in fasta format
-o (--ouput)	Output file names prefix [Default: genome]
-v (--vcf)	Input files in vcf format
-s (--step)	Size of the steps between windows [Default: 500]
-w (--window)	Width of the sliding windows [Default: 1000]
OPTIONS
die "$options\n" unless @ARGV;

my @vcf;
my @fasta;
my $output = 'genome';
my $slide = 500; my $window = 1000;

GetOptions(
	'v|vcf=s@{1,}' => \@vcf,
	'f|fasta=s@{1,}' => \@fasta,
	'o|ouput=s' => \$output,
	's|slide=i' => \$slide,
	'w|window=i' => \$window
);

## Populating list of SNPs
my %SNP;
my $tot_snps = 0;
while (my $vcf = shift@vcf){
	open IN, "<$vcf" or die "Can\'t open $vcf\n";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^#/){next;} # Skipping comments
		elsif ($line =~ /^(\S+)\s+(\d+).*HET=(\d);HOM=(\d)/){
			my $locus = $1;
			my $position = $2;
			my $het = $3; my $hom = $4; # HET or HOM
			if ($het != 0){$SNP{$locus}{$position} = 'het'; $tot_snps ++;}
			elsif ($hom != 0){$SNP{$locus}{$position} = 'hom'; $tot_snps ++;}
		}
	}
	close IN;
}

## Populating list of sequences
my %sequences;
my @seq_list;
my $seq;

while (my $fasta = shift @fasta) {
	open IN, "<$fasta" or die "Can\'t open $fasta\n";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+).*$/){$seq = $1; push (@seq_list, $seq);}
		else {$sequences{$seq} .= $line;}
	}
	close IN;
}

## Creating a "SNP" file for Circos
open SNP, ">$output.snp"; print SNP '#chr START END SNP_count'."\n";
my $mmax; my $mmin;

## Working on sequences
while (my $sequence = shift @seq_list){
	my $len = length($sequences{$sequence});
	my $terminus = $len - 1;
	my $max; my $min;
	print "Working on $sequence ... ";
	for(my $i = 0; $i <= $len - ($window-1); $i += $slide) {
		my $snp = 0;
		for ($i..($i+$window-1)){
			if (exists $SNP{$sequence}{$_}){$snp++;}
		}
		my $end = $i + $window - 1;
		print SNP "$sequence $i $end $snp\n";
		if (!defined $max){$max = $snp;}
		if (!defined $mmax){$mmax = $snp;}
		if (!defined $min){$min = $snp;}
		if (!defined $mmin){$mmin = $snp;}
		if ($snp > $max){$max = $snp;}
		if ($snp > $mmax){$mmax = $snp;}
		if ($snp < $min){$min = $snp;}
		if ($snp < $mmin){$mmin = $snp;}
	}
	print "Min/Max SNP per window = $min / $max\n";
}
print "Min SNP per window = $mmin \n";
print "Max SNP per window = $mmax \n";
print "Total number of SNPs = $tot_snps\n";
close SNP;
exit;
