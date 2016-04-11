#! /usr/bin/perl

use strict;
use warnings;

my $usage = 'perl parse_hits.pl hash query';

## Hash = blast hits
## Query = Helico_Protein.list

die $usage unless @ARGV;

open (HASH, $ARGV[0]);
my %hash = ();

#### Fill the hash with the hits.

	while (my $line = <HASH>) {
		chomp $line;
		
		if ($line =~ /^(H632_c\d{1,}\.\d{1,})/) {
			my $key = $1;
			my $value = $1;	## easier to look up if value = key
			$hash{"$key"} = "$value";
		}
	}

#### Work on the query file.

open (QUERY, $ARGV[1]);
open OUT, ">Coccomyxa.found";
open OUT2, ">Coccomyxa.notfound";

	while (my $line = <QUERY>){
		chomp $line;
		if ($line =~ /^(H632_c\d{1,}\.\d{1,})/){
			my $id = $1;
			if ($hash{$id} eq $id){
				print OUT "present in Coccomyxa\t$line\n";
			}else{
				print OUT2 "not found in Coccomyxa\t$line\n";
			}
		}
	}