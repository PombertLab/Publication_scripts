#!/usr/bin/env perl
## Pombert Lab, 2024

my $name = 'run_fastANI.pl';
my $version = 0.2;
my $updated = '2024-05-02';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage =<<"USAGE";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Runs fastANI in batch mode on specified FASTA files

REQS        fastANI - https://github.com/ParBLiSS/FastANI

COMMAND     ${name} \\
              -q *.fasta \\
              -r *.fasta \\
              -o ./fastANI \\
              -k 16 \\
              -t 8

OPTIONS:
-r (--ref)      FASTA file(s) to use as references
-q (--query)    FASTA file(s) to use as queries
-o (--outdir)   Output directory [Default: ./fastANI]
-k (--kmer)     Kmer to use [Default: 16]
-t (--threads)  CPU threads [Default: 8]
USAGE
die "\n$usage\n" unless @ARGV;

my @queries;
my @references;
my $outdir = './fastANI';
my $kmer = 16;
my $threads = 8;
GetOptions(
    'q|query=s@{1,}' => \@queries,
    'r|ref=s@{1,}' => \@references,
    'o|out=s' => \$outdir,
    'k|kmer=i' => \$kmer,
    't|threads=i' => \$threads
);

## Checking for fastANI
my $ani_check = `echo \$(command -v fastANI)`;
chomp $ani_check;
if ($ani_check eq ''){ 
    print STDERR "\n[E]: Cannot find fastANI. Please install fastANI in your \$PATH. Exiting..\n\n";
    exit;
}

## Creating output dir
unless (-d $outdir){
    mkdir ($outdir, 0755) or die "Can't create $outdir: $!\n";
}

## Creating reference list
my $ref_file = $outdir.'/references.list';
open RL, '>', $ref_file or die "Can't create $ref_file: $!\n";
foreach my $ref (@references){
    print RL $ref."\n";
}
close RL;

## Creating query list
my $que_file = $outdir.'/queries.list';
open QL, '>', $que_file or die "Can't create $que_file: $!\n";
foreach my $query (@queries){
    print QL $query."\n";
}
close QL;

## Running fastANI
my $outfile = $outdir.'/fastANI.results';
system (
    "fastANI \\
      --refList $ref_file \\
      --queryList $que_file \\
      --output $outfile \\
      --threads $threads \\
      --kmer $kmer \\
      --matrix
    "
);