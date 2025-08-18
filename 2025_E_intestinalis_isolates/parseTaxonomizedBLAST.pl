#!/usr/bin/env perl
## Pombert Lab, IIT 2019

my $name = 'parseTaxonomizedBLAST.pl';
my $version = '0.2b';
my $updated = '2024-06-02';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

################################################################################
# Usage definition
################################################################################

my $usage = <<"OPTIONS";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Parses the content of taxonomized BLAST searches

REQS        -outfmt '6 qseqid sseqid qstart qend pident length bitscore evalue staxids sscinames sskingdoms sblastnames'

USAGE       ${name} \\
              -b *.outfmt.6 \\
              -f *.fasta \\
              -n Streptococcus 'Streptococcus suis' 'Streptococcus sp.' \\
              -e 1e-25 \\
              -o output.fasta \\
              -v

OPTIONS:
-b (--blast)       BLAST input file(s)
-f (--fasta)       FASTA file(s)
-n (--name)        Names to be queried
-i (--inverse)     Invert selection; returns queries NOT matching specified names
-c (--column)      Which columns to query: sscinames, sskingdoms or sblastnames [Default: sscinames] 
-e (--evalue)      E-value cutoff for target organism(s) [Default: 1e-10]
-o (--output)      FASTA output file containing the desired sequences
-v (--verbose)     Verbose [Default: off]
--version          Show script version
OPTIONS

unless (@ARGV){
    print "\n$usage\n";
    exit(0);
};

################################################################################
### Command line options
################################################################################

my @blast;
my @fasta;
my @target;
my $column = 'sscinames';
my $evalue = 1e-10;
my $output;
my $inverse;
my $verbose;
my $sc_version;
GetOptions(
    'b|blast=s@{1,}' => \@blast,
    'f|fasta=s@{1,}' => \@fasta,
    'n|name=s@{1,}' => \@target,
    'c|column=s' => \$column,
    'e|evalue=s' => \$evalue,
    'o|output=s' => \$output,
    'I|inverse' => \$inverse,
    'v|verbose' => \$verbose,
    'version' => \$sc_version
);

unless ( ($column eq 'sscinames') || ($column eq 'sskingdoms') || ($column eq 'sblastnames') ){
    print "\nError. Column name $column is not recognised.\n";
    print "Please use either sscinames, sskingdoms or sblastnames\n\n";
    exit;
}

################################################################################
### Version
################################################################################

if ($sc_version){
    print "\n";
    print "Script:     $name\n";
    print "Version:    $version\n";
    print "Updated:    $updated\n\n";
    exit(0);
}

################################################################################
### Creating dababase of sequences from FASTA files
################################################################################

my %sequences;

for my $fasta (@fasta){
    open FASTA, "<", $fasta or die "Can't open file $fasta: $!\n";
    my $name;
    while (my $line = <FASTA>){
        chomp $line;
        if ($line =~ /^>(\S+)/){ $name = $1; }
        else { $sequences{$name} .= $line; }
    }
}

### Creating db of names to search for
my %scinames;
for my $names (@target){
    $scinames{$names} = $name;
}

################################################################################
### Parsing BLAST 'outfmt 6' file(s)
################################################################################

my %blasts;

while (my $blast = shift@blast){

    open BLAST, "<", $blast or die "Can't open file $blast: $!\n";

    while (my $line = <BLAST>){

        chomp $line;

        my @columns = split("\t", $line);
        my $query = $columns[0];
        my $bitscore = $columns[6];
        my $ev = $columns[7];
        
        ## Columns:
        ## [0] query, [1] target, [2] qstart, [3] qend, [4] pident, [5] length, 
        ## [6] bitscore, [7] evalue, [8] taxid, [9] sciname, [10] kingdom, [11] blastname
        
        unless ($ev <= $evalue){ next; }
        if (exists $blasts{$blast}{$query}){
            if ($bitscore > $blasts{$blast}{$query}[6]){ ## Checking for better hit(s) based on bitscores
                for (0..$#columns){
                    $blasts{$blast}{$query}[$_] = $columns[$_];
                }
            }
        }
        else {
            for (0..$#columns){
                $blasts{$blast}{$query}[$_] = $columns[$_];
            }
        }
    }
}

################################################################################
### Output file
################################################################################

open OUT, ">", $output or die "Can't create file $output: $!\n";
my @blasts = sort (keys %blasts); 

for my $blast (@blasts) {
    for my $query (sort (keys %{$blasts{$blast}})){
        
        ## Setting which column to parse
        my $regex;
        if ($column eq 'sscinames'){ $regex = $blasts{$blast}{$query}[9]; }
        elsif ($column eq 'sskingdoms'){ $regex = $blasts{$blast}{$query}[10]; }
        elsif ($column eq 'sblastnames'){ $regex = $blasts{$blast}{$query}[11]; }

        if ($verbose){ 
            print "Best hit for $query = $regex\n";
        }

        my @names = keys %scinames;
        my $flag = undef;
        for my $name (@names){
            if (!defined $flag){
                unless ($inverse){
                    if ($regex =~ /$name/i){
                        if ($verbose){
                            print "Match found for $name: ";
                            print "@{$blasts{$blast}{$query}}\n";
                        }
                        $flag = 'match';
                        fasta($query);
                    }
                }
                else {
                    if ($regex !~ /$name/i){
                        if ($verbose){
                            print "Match different from $name: ";
                            print "@{$blasts{$blast}{$query}}\n";
                        }
                        $flag = 'match';
                        fasta($query);
                    }
                }
            }
        }
    }
}

################################################################################
### Subroutine(s)
################################################################################

sub fasta {
    my $seq = $_[0];
    print OUT ">$seq\n";
    my @seq = unpack ("(A60)*", $sequences{$seq});
    while (my $tmp = shift@seq){
        print OUT "$tmp\n";
    }
} 