#!/usr/bin/perl

## Pombert Lab, IIT, 2014
## Runs Augutus on a set of fasta files
## Requires the augustus binary to be defined in the path

use strict;
use warnings;

my $usage = 'USAGE= run_augustus.pl HMMmodel *fastafiles';
die $usage unless @ARGV; 

my $hmm = shift@ARGV;	## Defining the hidden Markov model to be used with Augustus

while (my $file = shift@ARGV){
	if ($file =~ /.fsa$/){
		$file =~ s/.fsa$//;
		system 'augustus --gff3=on --species='."$hmm "."$file".'.fsa > '."$file".'.gff';
	}
	elsif ($file =~ /.fasta$/){
		$file =~ s/.fasta$//;
		system 'augustus --gff3=on --species='."$hmm "."$file".'.fsa > '."$file".'.gff';
	}
}
exit;