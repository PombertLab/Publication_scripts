#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2021

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;

my $name = 'vorocnn_average.pl';
my $version = 0.3;
my $update = '2021-10-22';

my $USAGE=<<"OPTIONS";

NAME    	$name
VERSION 	$version
UPDATE  	$update
SYNOPSIS	Calculate average VoroCNN scores per protein

EXAMPLE		${name} \\
		  -i RaptorX/*.scores AlphaFold/*.scores \\
		  -o vorocnn_averages.tsv

OPTIONS:
-i (--in)   	VoroCNN .scores file(s)
-o (--out)  	Output file name in .tsv format

OPTIONS

die "$USAGE" unless @ARGV;

my @VCNN; 
my $out = 'vorocnn_averages.tsv';

GetOptions(
    'i|in=s@{1,}' => \@VCNN,
    'o|out=s' => \$out
);

my %scores_db;
my %tools;

## Working on VoroCNN score files
while (my $file = shift@VCNN){
    open IN, "<", "$file";
    my ($base_name, $path) = fileparse($file);
    my $tool;
    if ($path =~ /alphafold/i){ $tool = 'alphafold'; $tools{$tool} = '';}
    elsif ($path =~ /raptorx/i){ $tool = 'raptorx'; $tools{$tool} = '';}
	elsif ($path =~ /refs/i){ $tool = 'RCSB references'; $tools{$tool} = '';}
    my ($locus_tag) = $base_name =~ /^(\w+)/;
    my $residue_count = 0;
    my $total_score = 0;

	while (my $line = <IN>) {
        chomp $line;
        if ($line =~ /^chain/){ next; }
        else{
            $residue_count++;
            my @columns = split("\t", $line);
            my $score = $columns[3];
            $total_score += $score;
        }
    }
    my $average = $total_score/$residue_count;
    $average = sprintf("%.3f", $average);
	$scores_db{$locus_tag}{$tool} = $average;

}

## Working on output TSV file for VoroCNN scores
open OUT, ">", "$out";
print OUT "Locus_tag";

my @tools = keys %tools; 
foreach my $tool (@tools){
	print OUT "\t$tool";
}
print OUT "\n";

foreach my $key (sort (keys %scores_db)){

	print OUT "$key";

	foreach my $tool (@tools){
		if (exists $scores_db{$key}{$tool}){
			print OUT "\t$scores_db{$key}{$tool}";
		}
		else { print OUT "\tN/A"; }
	}

	print OUT "\n";

}
close OUT;
