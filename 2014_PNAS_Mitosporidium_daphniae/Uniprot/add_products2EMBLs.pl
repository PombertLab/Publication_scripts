#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'add_products2EMBLs.pl products.list *.embl';
die $usage unless@ARGV;

### Filling the hash....
my $hh= shift@ARGV;
open HASH, "<$hh";
my %hash = ();

while (my $dbkey = <HASH>){
	chomp $dbkey;
	if ($dbkey =~ /^(\S+)\t(.*)$/){
		my $prot = $1;
		my $prod = $2;
		$hash{$prot}=$prod;
	}
}

### Working on the EMBL files...
while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.products";
	my $locus_tag = undef;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /FT\s+\/locus_tag=\"(DI09_\d+p\d+)\"/){
			$locus_tag=$1;
			print OUT "$line\n";
		}
		elsif ($line =~ /FT\s+CDS/){
			print OUT "$line\n";
			if (exists $hash{$locus_tag}){
				print OUT "FT                   /product=\"$hash{$locus_tag}\"\n";
			}
			else{
				print OUT "FT                   /product=\"hypothetical protein\"\n";;
			}
		}
		else{print OUT "$line\n";}
	}
}