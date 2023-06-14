#!/usr/bin/env perl
## Pombert Lab, Illinois Tech, 2023
my $name = 'parse_CVF_by_cc.pl'; 
my $version = '0.1';
my $updated = '2023-05-04';

use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(GetOptions);

my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Parses VCF files by start/end locations
		provided in a tab-delimited list

NOTE		cc => chromosome core; sub => (sub)telomeres

COMMAND		${name} \\
		  -l cc_list.txt \\
		  -v *.vcf \\
		  -o outdir

OPTIONS:
-l (--list)	Tab-delimited list: loci, start + end
-v (--vcf)	VCF files to parse
-o (--out)	Output directory [Default: ./]
USAGE
die "\n$usage\n" unless @ARGV;

my $cc_list;
my @vcfs;
my $outdir = './';
GetOptions(
	'l|list=s' => \$cc_list,
	'v|vcf=s@{1,}' => \@vcfs,
	'o|out=s' => \$outdir
);

## outdir
unless (-d $outdir){
	mkdir ($outdir, 0755) or die $!;
}

## CC
my %cc;
open CC, '<', $cc_list or die $!;

while (my $line = <CC>){
	
	chomp $line;
	
	unless ($line =~ /^#/){

		my @cols = split("\t", $line);
		my $locus = $cols[0];
		my $cc_start = $cols[1];
		my $cc_end= $cols[2];

		$cc{$locus}{'start'} = $cc_start;
		$cc{$locus}{'end'} = $cc_end;

	}
}

## VCFs
while (my $vcf = shift @vcfs){

	open VCF, '<', $vcf or die $!;

	my ($basename,$path) = fileparse($vcf);
	my ($prefix) = $basename =~ /^(.*)\.vcf$/;
	my $cc_outfile = $outdir.'/'.$prefix.'.cc.vcf';
	my $sub_outfile = $outdir.'/'.$prefix.'.sub.vcf';

	open my $cc, '>', $cc_outfile or die $!;
	open my $sub, '>', $sub_outfile or die $!;

	while (my $line = <VCF>){

		chomp $line;

		if ($line =~ /^#/){
			print $cc "$line\n";
			print $sub "$line\n";
		}
		else {

			my @cols = split("\t", $line); 
			my $locus = $cols[0];
			my $position = $cols[1];

			if ( ($position >= $cc{$locus}{'start'}) and ($position <= $cc{$locus}{'end'}) ){
				print $cc "$line\n";
			}
			else {
				print $sub "$line\n";
			}

		}

	}

}