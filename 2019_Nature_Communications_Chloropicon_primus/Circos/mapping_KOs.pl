#!/usr/bin/perl
## Pombert Lab, 2018
my $name = 'mapping_KOs.pl';
my $version = 0.1;

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Creates Circos tracks for KEGG orthologs
EXAMPLE		mapping_KOs.pl -l genome.list -c red -ko 2.1_Transcription.ko

OPTIONS:
-l	List of locus_tags and their positions for Circos ## See genome.list
-c	Desired color for the Circos highlights
-ko	KEGG orthology (.ko) reference file for desired pathway
OPTIONS
die "$options\n" unless @ARGV;


my $ko;
my $list;
my $color;
GetOptions(
	'ko=s' => \$ko,
	'l=s' => \$list,
	'c=s' => \$color
);
my %kos;
open IN, "<$ko";
while (my $line = <IN>){
	chomp $line;
	$kos{$line} = $line;
}
open LIST, "<$list";
open OUT, ">$ko.circos.kegg";
while (my $line = <LIST>){
	chomp $line;
	if ($line =~ /^\S+\t(\S+)\t(\d+)\t(\d+)\tplus\t(K\d+)$/){
		if (exists $kos{$4}){
			print OUT "$1 $2 $3 fill_color=$color,stroke_color=$color\n";
	}
	}
} 