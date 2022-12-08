#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2022
my $name = 'DUF.pl';
my $version = '0.1';
my $updated = '2022-07-11';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Checks for proteins with DUF (Domain of Unknown Function) annotations

COMMAND		${name} \\
		  -g input.gff \\
		  -o output.tsv \\
		  -p prefix \\
		  -r DUF

OPTIONS:
-g (--gff)	InterProScan GFF file
-o (--out)	Output file
-p (--prefix)	Protein prefix to search for [Default: 'JA071_\\w+']
-r (--regex)	Pattern/regex to search for [Default: 'DUF\\d+']
USAGE
die "\n$usage\n" unless @ARGV;

my $gff;
my $out;
my $prefix = 'JA071_\w+';
my $regex = 'DUF\d+';
GetOptions(
	'g|gff=s' => \$gff,
	'o|out=s' => \$out,
	'p|prefix=s' => \$prefix,
	'r|regex=s' => \$regex
);

open GFF, "<", $gff or die "Can't open $gff: $!\n";
open OUT, ">", $out or die "Can't create $out: $!\n";

my %matches;
while (my $line = <GFF>){
	chomp $line;
	if ($line =~ /^($prefix).*($regex)/){
		push (@{$matches{$1}}, $2);
	}
}

for my $key (sort (keys %matches)){
	my @array = @{$matches{$key}};
	if (scalar @array == 1){
		print OUT "$key\t$array[0]\n";
	}
	else {
		print OUT "$key\t";
		for (0..$#array){
			print OUT "$array[$_], ";
		}
		print OUT "$array[-1]\n";
	}
}