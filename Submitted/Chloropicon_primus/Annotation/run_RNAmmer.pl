#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perk run_RNAmmer.pl kingdom files.fsa'; ## kingdom = arc, bac or euk

my $kingdom = shift@ARGV;

while (my $file = shift@ARGV){
	system "/opt/RNAmmer/rnammer -S $kingdom -m tsu,ssu,lsu -gff $file.gff2 -h $file.hmm -f $file.rRNAs < $file";
}