#!/usr/bin/env perl
## Pombert Lab, 2025

my $name = 'cafe_sig_tree.pl';
my $version = '0.1a';
my $updated = '2025-06-18';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);

###################################################################################################
### Command line options
###################################################################################################

my $usage = <<"USAGE";
NAME            ${name}
VERSION         ${version}
UPDATED         ${updated}
SYNOPSIS        Extract trees with significant p-values from CAFE5 *_asr.tre file and create lists
                of significant changes per nodes

EXAMPLE:        ${name} \\
                  -i Base_asr.tre \\
                  -t Base_change.tab \\
                  -i Base_count.tab \\
                  -r Base_branch_probabilities.tab \\
                  -p significant_trees

OPTIONS:
-i (--input)    CAFE5 Base_asr.tre file
-t (--tab)      CAFE5 Base_change.tab file
-u (--count)    CAFE5 Base_count.tab file
-r (--prob)     CAFE5 Base_branch_probabilities.tab file
-d (--dir)      Desired output directory [Default: SIGS]
-p (--prefix)   Desired output files prefix
-v (--version)  Show script version
USAGE

unless (@ARGV){
    print "\n$usage\n";
    exit(0);
};

my @commands = @ARGV;
my $cafe_tree;
my $change_tab;
my $count_tab;
my $prob_tab;
my $outdir = 'SIGS';
my $prefix;
my $sc_version;

GetOptions(
    'i|input=s' => \$cafe_tree,
    't|tab=s' => \$change_tab,
    'u|count=s' => \$count_tab,
    'r|prob=s' => \$prob_tab,
    'o|outdir=s' => \$outdir,
    'p|prefix=s' => \$prefix,
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
### Outdir
#########################################################################

unless (-d $outdir) {
    make_path($outdir,{mode => 0755}) or die "Can't create $outdir: $!\n";
}

#########################################################################
### Create change/count/probabilities databases

open CH, '<', $change_tab or die "Can't read $change_tab: $!\n";
open CT, '<', $count_tab or die "Can't read $count_tab: $!\n";
open PB, '<', $prob_tab or die "Can't read $prob_tab: $!\n";

my %change_db;
my %count_db;
my %prob_db;

my $nodes_number;

# changes
while (my $line = <CH>){

    chomp $line;

    if ($line =~ /^FamilyID/){
        my @data = split("\t", $line);
        $nodes_number = scalar(@data) - 1;
    }
    else{
        my @data = split("\t", $line);
        for my $num (1..$#data){
            $change_db{$data[0]}{$num} = $data[$num];
        }

    }


}

# counts
while (my $line = <CT>){

    chomp $line;

    if ($line =~ /^FamilyID/){
        next;
    }
    else{
        my @data = split("\t", $line);
        for my $num (1..$#data){
            $count_db{$data[0]}{$num} = $data[$num];
        }

    }

}

# probabilities
while (my $line = <PB>){

    chomp $line;

    if ($line =~ /^FamilyID/){
        next;
    }
    else{
        my @data = split("\t", $line);
        for my $num (1..$#data){
            $prob_db{$data[0]}{$num} = $data[$num];
        }

    }

}

close CH;
close CT;
close PB;

#########################################################################
### Parse significant trees

for my $digit (1..$nodes_number){

    my $tree_file = $outdir.'/'.$prefix.'.c'.$digit.'.tre';
    my $list_file = $outdir.'/'.$prefix.'.c'.$digit.'.list';

    open IN, '<', $cafe_tree or die "Can't read $cafe_tree: $!\n";
    open TREE, '>', $tree_file or die "Can't create $tree_file: $!\n";
    open LIST, '>', $list_file or die "Can't create $list_file: $!\n";

    print LIST '# Orthogroup'."\t"."Change"."\t".'Count'."\t".'P-value'."\n";

    while (my $line = <IN>){
        chomp $line;
        if ($line !~ /^  TREE /){
            print TREE $line."\n";
        }
        else{
            if ($line =~ /^  TREE (\S+).*\<($digit)\>\*/){

                my $og = $1;
                print LIST $og."\t".$change_db{$og}{$digit};
                print LIST "\t".$count_db{$og}{$digit};
                print LIST "\t".$prob_db{$og}{$digit}."\n";
                print TREE $line."\n";

            }
        }
    }

    close IN;
    close TREE;
    close LIST;

}
