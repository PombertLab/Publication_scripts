#!/usr/bin/env perl
## Pombert Lab, 2025

my $name = 'expand_og_list.pl';
my $version = '0.1';
my $updated = '2025-06-19';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

###################################################################################################
### Command line options
###################################################################################################

my $usage = <<"USAGE";
NAME            ${name}
VERSION         ${version}
UPDATED         ${updated}
SYNOPSIS        Expand OrthoFinder orthogroups and add their function from SYNY annotation files

EXAMPLE:        ${name} \\
                  -i Orthogroups.txt \\
                  -a *.lists \\
                  -o expanded_ogs.txt

OPTIONS:
-i (--input)    OrthoFinder Orthogroups.txt file
-a (--annot)    SYNY annotation (.list) file(s)
-o (--out)      Desired output file
-v (--version)  Show script version
USAGE

unless (@ARGV){
    print "\n$usage\n";
    exit(0);
};

my @commands = @ARGV;
my $ogs;
my @annots;
my $outfile;
my $sc_version;

GetOptions(
    'i|input=s' => \$ogs,
    'a|annot=s@{1,}' => \@annots,
    'o|out=s' => \$outfile,
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
### Creating annotations database
#########################################################################

my %annotations;

for my $annot (@annots){

    open AN, '<', $annot or die "Can't open $annot: $!\n";

    while (my $line = <AN>){

        chomp $line;

        if ($line =~ /^$/){
            next;
        }
        else{
            my @data = split("\t", $line);
            my $locus = $data[0];
            my $annotation = $data[-1];
            $annotations{$locus} = $annotation;
        }

    }

    close AN;

}

#########################################################################
### Creating expanded orthogroups file
#########################################################################

open IN, '<', $ogs or die "Can't read $ogs: $!\n";
open OUT, '>', $outfile or die "Can't create $outfile: $!\n";

while (my $line = <IN>){

    chomp $line;

        if ($line =~ /^$/){
        next;
    }
    else{

        my @data = split(' ', $line);
        my $og = $data[0];
        $og =~ s/:$//;

        print OUT '### '.$og."\n";
        for my $num (1..$#data){
            print OUT '  '.$data[$num]."\t".$annotations{$data[$num]}."\n";
        }
        print OUT "\n";
    
    }

}