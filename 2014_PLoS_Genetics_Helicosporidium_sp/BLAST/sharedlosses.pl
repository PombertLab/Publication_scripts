#! /usr/bin/perl

use strict;
use warnings;

my $usage = 'perl shared.pl hash query';

## Hash = blast hits
## Query = Helico_Protein.list

die $usage unless @ARGV;

open (HASH, $ARGV[0]);
my %hash = ();

#### Fill the hash with the hits.

	while (my $line = <HASH>) {
		chomp $line;
		
		if ($line =~ /^not found in Helico\t(\S+)/) {
			my $key = $1;
			my $value = $1;	## easier to look up if value = key
			$hash{"$key"} = "$value";
		}
	}

#### Work on the query file.

open (QUERY, $ARGV[1]);
open OUT, ">missingfromboth.txt";
open OUT2, ">unique.txt";

	while (my $line = <QUERY>){
		chomp $line;
		if ($line =~ /^not found in Helico\t(\S+)/){
			my $id = $1;
			if ($hash{$id} eq $id){
				print OUT "$line\n";
			}else{
				print OUT2 "$line\n";
			}
		}
	}