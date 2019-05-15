#!/usr/bin/perl
## Pombert Lab, 2017
my $version = '0.1';
my $name = 'getProducts.pl';


## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Generates a product list from the NCBI protein fasta files (*.faa)
USAGE		getProducts.pl *.faa
OPTIONS
die "$options\n" unless @ARGV;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.faa$//;
	open OUT, ">$file.products";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+)\s+(.*)\s+\[/){
			my $locus = $1;
			my $product = $2;
			print OUT "$locus\t$product\n";
		}
	}
}