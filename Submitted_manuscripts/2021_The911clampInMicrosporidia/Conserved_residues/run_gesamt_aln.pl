#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;

my $name = 'run_gesamt_aln.pl';
my $version = '0.1b';
my $update = '2022-03-14';

my $usage=<<"USAGE";

NAME            $name
VERSION         $version
UPDATE          $update
SYNOPSIS        Pairwise 3D structure alignment with GESAMT

EXAMPLE         ./${name} -q *.pdb -t *.pdb -a

OPTIONS
-q (--query)    Proteins to use as query (fixed)
-t (--target)   Proteins to use as target (moving)
-o (--output)   Prefix for output files
-a (--align)    To generate alignment files in fasta and txt format
-p (--pdb)      To generate alignment files plus superimposed PDB

USAGE

my @query;
my @target;
my $align;
my $pdb;
my $out = 'out';

GetOptions(

    'q|query=s@{1,}' => \@query,
    't|target=s@{1,}' => \@target,
    'a|align' => \$align,
    'p|pdb' => \$pdb,
    'o|output=s' => \$out
);

# Query files
while (my $file1 = shift@query){

    my $file1_name = fileparse($file1);
    
    # For each query to run against every target
    for (@target){
        my $file2 = $_;
        my $file2_name = fileparse($file2);

        # Generating alignment and txt files with GESAMT
        if ($align){
            print "Performing alignment only on $file1_name and $file2_name...\n";
            system("gesamt \\
                    $file1 \\
                    $file2 \\
                    -a ${out}.fasta \\
                    > ${out}.txt"
            );
        }
        # Generating alignment, txt and pdb of superimposed structures with GESAMT
        elsif ($pdb){
            print "Performing alignment and generating PDB on $file1_name and $file2_name...\n";
            system("gesamt \\
                    $file1 \\
                    $file2 \\
                    -o ${out}.pdb \\
                    -a ${out}.fasta \\
                    > ${out}.txt"
            );
        }
    }
}