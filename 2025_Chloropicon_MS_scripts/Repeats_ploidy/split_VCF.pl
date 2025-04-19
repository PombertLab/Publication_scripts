#!/usr/bin/env perl
# Pombert lab 2019

my $version = '0.3';
my $name = 'split_VCF.pl';
my $updated = '2025-01-14';

###################################################################################################
# Pragmas / required modules
###################################################################################################

use strict;
use warnings;
use PerlIO::gzip;
use File::Basename;
use Getopt::Long qw(GetOptions);

###################################################################################################
# README / Command line options
###################################################################################################

my $options = <<"OPTIONS";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Splits VCF file(s) per contig/chromosome. Can also exclude masked positions (N) 
            in the reference (e.g. by tools such as RepeatMasker).

NOTE        Assumes that contigs/chromosomes have unique names (e.g. accession numbers). VCF
            files with conflicting sequence names should be handled separately.

COMMAND     ${name} \\
              --vcf *.vcf \\
              --exclude \\
              --outdir ./SPLIT

OPTIONS:
-v (--vcf)      VCF file(s) to parse (GZIP files are supported)
-o (--outdir)   Output directory [Default: SPLIT]
-e (--exclude)  Exclude masked (N) positions
-n (--nohead)   Exclude VCF header from parsed output files
-version        Show script version
OPTIONS

unless (@ARGV){
	print "\n$options\n";
	exit(0);
};

my @VCFs;
my $outdir = 'SPLIT';
my $exclude;
my $nohead;
my $scversion;
GetOptions(
    'v|vcf=s@{1,}' => \@VCFs,
    'o|outdir=s' => \$outdir,
    'e|exclude' => \$exclude,
    'n|nohead' => \$nohead,
    'version' => \$scversion
);

###################################################################################################
# Version
###################################################################################################

if ($scversion){
    print("\n");
    print('Script:  '.$name."\n");
    print('Version: '.$version."\n");
    print('Updated: '.$updated."\n\n");
    exit(0);
}

###################################################################################################
# Output dir
###################################################################################################

unless (-d $outdir){
    mkdir ($outdir, 0755) or die $!;
}

###################################################################################################
# Output VCF files
###################################################################################################

while (my $vcf = shift@VCFs){

    my $file_name = basename($vcf);
    my @file_data = split('\.',$file_name);

    my $diamond = "<";
    if ($file_data[-1] eq 'gz'){ 
        $diamond = "<:gzip";
    }

    open VCF, $diamond, $vcf or die "Can't read $vcf: $!\n";

    my @header;
    my %variants;

    my $previous_contig = undef;
    my $filename;

    while (my $line = <VCF>){

        chomp $line;

        ## Data header
        if ($line =~ /^#/){
            push(@header, $line);
        }

        ## Parsing VCF data
        else {

            my @columns = split("\t", $line);
            my $contig = $columns[0];
            my $ref_nt = $columns[3];

            if ($exclude){
                if ($ref_nt eq 'N'){
                    next;
                }
                else {
                    push(@{$variants{$contig}}, $line);
                }
            }
            else {
                push(@{$variants{$contig}}, $line)
            }

        }

    }

    ## Splitting data per contig/chromosome
    for my $contig (sort(keys %variants)){

        # Progess report
        print 'Parsing: '.$vcf.' : '.$contig."\n";

        my $filename = $outdir.'/'.$contig.'.vcf';
        open my $fh, '>', $filename or die "Can't create $contig: $!\n";

        unless ($nohead){
            for my $line (@header){
                print $fh $line."\n";
            }
        }

        for my $line (@{$variants{$contig}}){
            print $fh $line."\n";
        }

        close $fh;

    }

}
