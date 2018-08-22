#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = 'USAGE = mapping_KOs.pl -l genome.list -c red -ko 2.1_Transcription.ko';
die "\n$usage\n\n" unless @ARGV;

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