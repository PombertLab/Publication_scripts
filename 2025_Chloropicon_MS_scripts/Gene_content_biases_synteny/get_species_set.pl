#!/usr/bin/env perl
# Pombert Lab, 2023

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Copy;

my $name = 'get_species_set.pl';
my $version = '0.2';
my $updated = '2023-11-19';

my $usage = << "USAGE";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS

COMMAND     ${name} \\
              -l *.txt \\
              -x RCC1871 CCMP1998 RCC2335 \\
              -f ./OrthoFinder/Results/Results_Nov16/Orthogroup_Sequences \\
              -d ./OrthoFinder/Results/Results_Nov16/Orthogroups/Orthogroups.tsv \\
              -a ../SYNY/ALL/ANNOTATIONS/*.annotations \\
              -outdir C_primus \\
              -prefix C_primus_og

OPTIONS:
-l (--list)     List of orthogroups generated with parse_OGs.pl
-f (--fasta)    Folder containing the orthogroup FASTA sequences
-d (--desc)     Orthogroups content description (Orthogroups.tsv)
-a (--annots)   List of products and their annotations generated with SYNY
-x (--exclude)  Name of strains to exclude
-o (--outdir)   Output directory [Default: ./]
-p (--prefix)   Output file prefix [Default: og_set]
USAGE
die "\n$usage\n" unless @ARGV;

my @lists;
my $og_fasta;
my @og_desc;
my @annots;
my $outdir = './';
my $prefix = 'og_set';
my @excluded;
GetOptions(
    'l|list=s@{1,}' => \@lists,
    'f|fasta=s' => \$og_fasta,
    'd|desc=s@{1,}' => \@og_desc,
    'o|outdir=s' => \$outdir,
    'a|annot=s@{1,}' => \@annots,
    'p|prefix=s' => \$prefix,
    'x|exclude=s@{1,}' => \@excluded
);

## Checking outdir
unless (-d $outdir){
    mkdir ($outdir, 0755) or die "Can't create $outdir: $!\n";
}

my $fasta_dir = $outdir.'/'.'FASTA';
unless (-d $fasta_dir){
    mkdir ($fasta_dir, 0755) or die "Can't create $fasta_dir: $!\n";
}

my $outfile = $outdir.'/'.$prefix.'.txt';
my $prod_file = $outdir.'/'.$prefix.'.products.txt';

## Parsing annotations
my %annots;

while (my $file = shift @annots){
    open ANN, '<', $file or die "Can't open $file: $!\n";
    my ($basename) = fileparse($file);
    $basename =~ s/\.\w+$//;
    while (my $line = <ANN>){
        chomp $line;
        my ($locus, $product) = split("\t", $line);
        $annots{$locus} = $product;
    }
}

## Parsing Orthogroups.tsv
my %og_desc;
while (my $og_desc = shift @og_desc){
    open OGG, '<', $og_desc or die "Can't open $og_desc: \n";
    while (my $line = <OGG>){
        chomp $line;
        unless ($line =~ /^Orthogroup/){
            $line =~ s/\r$//; ## In case of MSDOS format
            my @data = split("\t", $line);
            my $og = shift @data;
            foreach my $entry (@data){
                my @items = split(", ", $entry);
                for (@items){
                    $og_desc{$og}{$_} = 1;
                }
            }
        }
    }
}

## Parsing lists
my %species;
my %orthogroups;
my %excluded;

for my $ex (@excluded){
    $excluded{$ex} = 1;
}

while (my $list = shift @lists){

    open LIST, '<', $list or die "Can't open $list: $!\n";
    my ($basename) = fileparse($list);
    $basename =~ s/\.(\w+)$//;

    while (my $line = <LIST>){
        chomp $line;
        $species{$basename}{$line} = 1;
        $orthogroups{$line} = 1;
    }

}

## Printing sets that are not exluded...

open OUT, '>', $outfile or die "Can't create $outfile: $!\n";
open PRD, '>', $prod_file or die "Can't create $prod_file: $!\n";

foreach my $og (sort (keys %orthogroups)){

    my $exclude_counter = 0;

    foreach my $ex (keys %excluded){
        if (exists $species{$ex}{$og}){
            $exclude_counter = 1;
        }
    }

    if ($exclude_counter == 0){

        print OUT $og."\n";
        print PRD '### '.$og."\n";

        foreach my $key (keys %{$og_desc{$og}}){
            my $product = $annots{$key};
            print PRD '  '.$key."\t".$product."\n";
        }
        print PRD "\n";

        ## Copy file to subdir
        my $source_file = $og_fasta.'/'.$og.'.fa';
        my $dest_file = $fasta_dir.'/'.$og.'.fa';
        copy($source_file, $dest_file);

    }

}



