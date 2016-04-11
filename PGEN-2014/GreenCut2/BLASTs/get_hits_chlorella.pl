#! /usr/bin/perl

use strict;
use warnings;

## Hash = chlamy.3sources.filtered --> i.e. filtered BLAST output
## Query = the suplplemental table jbc.M111.233734-2 from the GreenCut2 paper as text file

my $usage = 'perl get_hits.pl hash query';
die $usage unless @ARGV;
open (HASH, $ARGV[0]);
my %hash = ();

#### Fill the hash with the hits.
while (my $line = <HASH>) {
	chomp $line;
	if ($line =~ /^(\d{1,})\t/) {
		my $key = $1;
		my $value = $1;	## easier to look up if value = key
		$hash{"$key"} = "$value";
	}
}
#### Work on the query file.
open (QUERY, $ARGV[1]);
open OUT, ">Chlorella.found";
open OUT2, ">Chlorella.notfound";
while (my $line = <QUERY>){
	chomp $line;
	if ($line =~ /^(\d{1,})\t/){
		my $id = $1;
		if ($hash{$id} eq $id){
			print OUT "present in Chlorella\t$line\n";
		}else{
			print OUT2 "not found in Chlorella\t$line\n";
		}
	}
}