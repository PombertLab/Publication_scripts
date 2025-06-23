#!/usr/bin/env perl
## Pombert Lab, 2025

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $name = 'og_to_cafe.pl';
my $version = '0.1';
my $updated = '2025-06-17';

###################################################################################################
### Command line options
###################################################################################################

my $usage = <<"USAGE";
NAME            ${name}
VERSION         ${version}
UPDATED         ${updated}
SYNOPSIS        Converts OrthoFinder's Orthogroups.GeneCount.tsv file to format compatible for
                CAFE5 analyses

EXAMPLE:        ${name} -o Orthogroups.GeneCount.tsv -c cafe.tsv

OPTIONS:
-o (--og)       Orthogroups.GeneCount.tsv input file from OrthoFinder
-c (--cafe)     Desired output file (in CAFE5 format)
-a (--anno)     Tab-delimited orthogroup/annotation (Optional)
-v (--version)  Show script version
USAGE

unless (@ARGV){
    print "\n$usage\n";
    exit(0);
};

my @commands = @ARGV;
my $og_tsv;
my $cafe_tsv;
my $annots;
my $sc_version;

GetOptions(
    'o|og=s' => \$og_tsv,
    'c|cafe=s' => \$cafe_tsv,
    'a|anno=s' => \$annots,
    'v|version' => \$sc_version
);

#########################################################################
### Version
#########################################################################

if ($sc_version){
    print "\n";
    print "Script:     $name\n";
    print "Version:    $version\n";
    print "Updated:    $updated\n\n";
    exit(0);
}

#########################################################################
### Load orthogroup annotations
#########################################################################

my %og_annots;

if ($annots){

    open AN, '<', $annots or die "Can't read $annots: $!\n";
    ## Data stucture => Orthogroup name\tOrthogroup annotation\n

    while (my $line = <AN>){

        chomp $line;

        if ($line =~ /^#/){
            next; ## Skip comments 
        }
        else{
            my @data = split("\t", $line);
            $og_annots{$data[0]} = $data[1];
        }

    }

}

#########################################################################
### Converting OG to CAFE5
#########################################################################

open OG, '<', $og_tsv or die "Can't read $og_tsv: $!\n";
open CF, '>', $cafe_tsv or die "Can't create $cafe_tsv: $!\n";

while (my $line = <OG>){

    chomp $line;

    my @data = split("\t", $line);
    pop @data; ## Removing the total column

    if ($data[0] eq 'Orthogroup'){
        print CF "Desc\t";
    }
    else{
        if (exists $og_annots{$data[0]}){
            print CF $og_annots{$data[0]}."\t";
        }
        else{
            print CF '(null)'."\t";
        }
    }

    for (0..$#data-1){
        print CF $data[$_]."\t";
    }
    print CF $data[-1]."\n";

}