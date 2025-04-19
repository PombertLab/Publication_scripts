#!/usr/bin/env perl
## Pombert Lab, 2023

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $name = 'parse_OGs.pl';
my $version = 0.1;
my $updated = '2023-11-17';

my $usage = << "USAGE";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Creates lists of orthogroups found per species. Useful
            to generate Venn diagrams with 
            https://bioinformatics.psb.ugent.be/webtools/Venn/

COMMAND     ${name} \\
              -og *.tsv \\
              -out FOLDER/

OPTIONS:
-og      Orthogroups.tsv + Orthogroups_UnassignedGenes.tsv from OrthoFinder
-out     Output folder [Default: ./]
USAGE
die "\n$usage\n" unless @ARGV;

my @ogs;
my $outdir = './';
GetOptions(
    'og=s@{1,}' => \@ogs,
    'out=s' => \$outdir
);

## Checking output directory
unless (-d $outdir){
    mkdir ($outdir, 0755) or die "Can't create $outdir: $!\n";
}

## Parsing the OrthoFinder files

my %organism;
my %position;

for my $og (@ogs){

    open OG, '<', $og or die $!;

    while (my $line = <OG>){

        chomp $line;

        my @data = split("\t", $line);

        if ($line =~ /^Orthogroup/){
            for (1..$#data){
                $position{$_} = $data[$_];
                $organism{$data[$_]} = $data[$_];
            }
        }

        else {
            my $og_group = $data[0];
            for (1..$#data){
                my $species = $position{$_};
                if ($data[$_] ne ''){
                    no strict;
                    push (@{$organism{$species}}, $og_group);
                }
            }
        }
    }
}

## Printing data for each organism
foreach my $key (keys %organism){
    open FH, '>', "$outdir/$key.txt" or die "Can't create $outdir/$key.txt: $!\n";
    no strict;
    while (my $og = shift @{$organism{$key}}){
        print FH $og."\n";
    }
}