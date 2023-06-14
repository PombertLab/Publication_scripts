#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(GetOptions);

die "\nUSAGE = $0 -v *.vcf -o outdir\n" unless @ARGV;

my @vcf;
my $outdir = './';
GetOptions(
	'v|vcf=s@{1,}' => \@vcf,
	'o|outdir=s' => \$outdir
);

unless (-d $outdir){
	mkdir ($outdir, 0755) or die $!;
}

while (my $vcf = shift@vcf){
	my ($basename, $path) = fileparse($vcf);
	my ($query,$ref) = $basename =~ /^(\w+).*.(?:R1|SE).fastq.gz.(\w+).*.vcf$/;
	my $outfile = $query.'.vs.'.$ref.'.vcf';
	# print "$outfile\n";
	system "cp $vcf $outdir/$outfile";
}