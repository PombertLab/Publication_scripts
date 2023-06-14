#!/usr/bin/perl
## Carol, Pombert Lab, IIT 2022

use strict;
use warnings;

my $name = 'revcomp.pl';
my $version = '0.1';
my $update = '2022-12-19';

my $usage = <<"USAGE";

NAME        $name
VERSION     $version
UPDATE      $update
SYNOPSIS    Reverse complement DNA sequences in fasta file format.

EXAMPLE     ./$name *.fasta

USAGE
die "$usage" unless @ARGV;

while (my $file = shift@ARGV){

    open IN, "<", "$file";
    $file =~ s/.fasta/_reverse.fasta/;
    $file =~ s/.fsa/_reverse.fsa/;
    open OUT, ">", "$file";

    my %chrom;
    my $name;

    while (my $line = <IN>){
        chomp $line;

        if ($line =~ /^>(\w+)/ ){
            $name = $1;
        }
        else {
            $chrom{$name} .= $line;
        }
    }

    for (sort (keys %chrom)){

        print OUT ">$_\n";
        my $seq = $chrom{$_};
        my $reverse = reverse($seq);                                            # Reverse orientation
        $reverse =~ tr/AaCcGgTtRrYyKkMmBbVvDdHh/TtGgCcAaYyRrMmKkVvBbHhDd/;      # Complement of seq
        my @fna = unpack("(A60)*", $reverse);

        while (my $nucl = shift@fna){
            print OUT "$nucl\n";
        }
    }
}