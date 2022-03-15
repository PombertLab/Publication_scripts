#!/usr/bin/perl
# Pombert Lab, Illinois Institute of Technology, 2022

use strict;
use warnings;

my $name = '3dalign_heatmap.pl';
my $version = '0.2a';
my $update = '2022-03-07';

my $usage = <<"USAGE";

NAME        ${name}
VERSION     ${version}
UPDATE      ${update}
SYNOPSIS    Generates per-residue similarity score in .csv file and alignment .stats file from GESAMT alignment

EXAMPLE ./${name} *.txt

USAGE

die "$usage" unless @ARGV;

# Assigning similarity score to GESAMT similarity level symbols
my %symbols = (
    '**' => 5,
    '++' => 4,
    '==' => 3,
    '--' => 2,
    '::' => 1,
    '..' => 0
);

my %hm_data;

while (my $file = shift@ARGV){

    open IN, "<", "$file";
    $file =~ s/.txt//;
    open OUT, ">", "${file}_aln.csv";
    open STAT, ">", "${file}_aln.stats";

    my @length;
    my $length;
    my $aligned_residue;
    my $similar_count = 0;
    my $aa_count = 0;

    # Printing header in csv file
    print OUT "AA,Identity\n";

    while (my $line = <IN>){
        chomp $line;
        if ($line =~ /^\|\s+FIXED/){
            next;
        }
        # Grabbing number of aligned residues
        elsif($line =~ /^\s+Aligned\s+residues\s+:\s+(\d+)/){
            $aligned_residue = $1;
        }

        # Grabbing identity and similarity info from alignment
        elsif ($line =~ /^\|(.*?)\|(.*?)\|(.*?)\|/){
            my $fixed = $1;
            my $dist = $2;
            my $mov = $3;
            
            ## Working on fixed
            if ($fixed =~ /(\w+)\s+(\d+)\s+$/){
                my $aa = $1;
                my $pos = $2;

                # Grabbing length of protein
                push(@length, $pos);
                $length = scalar(@length);

                # Building score list based on similarity level
                if ($dist =~ /^\s+\<(..)/){
                    my $symb = $1;
                    my $number = $symbols{$symb};
                    print OUT "${aa}_$pos,$number\n";

                    # Counting identical amino acids in alignment
                    if ($symb =~ /\*\*/){
                        $aa_count++;
                    }

                    # Counting similar amino acids in alignment
                    if($symb ne '..'){
                        $similar_count++;
                    }
                }
                # Assigning score penalty for gaps
                else {
                    my $number = -1;
                    print OUT "${aa}_$pos,$number\n";
                }
            }

        }
    }

    # Calculating amino acid identity in %
    # my $stats_identity = sprintf("%.2f", ($seq_identity * 100));
    my $stats_identity = sprintf("%.2f", ($aa_count/$length * 100));
    my $stats_similarity = sprintf("%.2f", ($similar_count/$aligned_residue *100));

    # Populating heatmap db
    $hm_data{$file} = $stats_similarity;

    # Printing to the output stats file
    print STAT "Query protein length = $length\n";
    print STAT "Aligned residues = $aligned_residue\n";
    print STAT "Number of identical residues = $aa_count\n";
    print STAT "Number of similar residues = $similar_count\n";
    print STAT "Pairwise identity (%) = $stats_identity\n";
    print STAT "Pairwise similarity (%) = $stats_similarity\n";

}