#!/usr/bin/perl
## Pombert Lab, Illinois Tech 2022
my $name = 'split_GBK.pl';
my $version = '0.1a';
my $updated = '2022-10-20';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);

my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Splits a GBK file into multiple GBK files, one per contig/locus

COMMAND		${name} -g file.gbk -o output_dir -e gbk

OPTIONS:
-g (--gbk)	Input GenBank file (.gb|gbff|gbk) file to split
-o (--out)	Output directory [Default: ./]
-e (--ext)	Desired file extension [Default: gbk]
USAGE
die "\n$usage\n" unless @ARGV;

my $gbk;
my $outdir = './';
my $ext = 'gbk';
GetOptions(
	'g|gbk=s' => \$gbk,
	'o|out=s' => \$outdir,
	'e|ext=s' => \$ext
);

### Checking output directory
unless (-d $outdir){
	make_path($outdir, { verbose => 0, mode => 0755}) or die "Can't create $outdir: $!\n";
}

### Creating a single string per GBK, then splitting per // 
open IN, "<", $gbk or die "Can't open $gbk: $!\n";
while (my $line = <IN>){
	$gbk .= $line;
}
close IN;

my @contigs = split(/\/\//, $gbk);

### Writing separate GBK files
foreach my $contig (@contigs){
	if ($contig  =~ /LOCUS       (\S+)/){
		my ($locus) = $contig =~ /LOCUS       (\S+)/;
		my $outfile = $outdir.'/'.$locus.'.'.$ext;
		open OUT, ">", $outfile or die "Can't create $outfile: $!\n";
		print OUT "$contig";
		print OUT "//\n";
	}
}
