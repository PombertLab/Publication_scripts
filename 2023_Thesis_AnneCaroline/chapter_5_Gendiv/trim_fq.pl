#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2023

my $name = 'trim_fq.pl';
my $version = '0.1a';
my $updated = '2023-05-22';

use strict;
use warnings;
use File::Basename;
use PerlIO::gzip;
use Getopt::Long qw(GetOptions);

my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UDPATED		${updated}
SYNOPSIS	Trim FASTQ reads to desired length (keeps first X bases)

COMMAND		${name} \\
		  -f *.fastq \\
		  -l 100 \\
		  -n 1 \\
		  -o DATA

OPTIONS:
-f (--fastq)	FASTQ file(s) to trim
-l (--len)		First X bases to keep [Default: 100]
-n (--num)		Number of reads to keep (in millions) [Default: 1]
-o (--out)		Output directory [Default: ./]
USAGE
die "\n$usage\n" unless @ARGV;

my @fastq;
my $len = 100;
my $num = 1;
my $outdir = './';
GetOptions(
	'f|fastq=s@{1,}' => \@fastq,
	'l|len=i' => \$len,
	'n|num=i' => \$num,
	'o|out=s' => \$outdir
);

unless (-d $outdir){
	mkdir ($outdir, 0755) or die $!;
}


while (my $fastq = shift @fastq){

	my $counter = 0;
	
	## Autodetecting if file is gzipped from the file extension

	my $FQ;
	my $format;

	if ($fastq =~ /.gz$/){
		open $FQ, '<:gzip', $fastq or die $!;
		$format = 'gzip';
	}
	else {
		open $FQ, '<', $fastq or die $!;
		$format = 'fastq';
	}

	my ($basename) = fileparse($fastq);
	$basename =~ s/\.fastq.gz$//;
	$basename =~ s/\.fastq$//;

	open my $fh, '>', "$outdir/$basename.trimmed.${num}M.fastq" or die $!;

	my $total = $num * 1000000 * 4; 

	while (my $line = <$FQ>){
		
		chomp $line;
		
		$counter++;
		my $modulo = $counter % 4;

		# Checking if line is sequence or quality score and trim
		if (($modulo == 2) or ($modulo == 0)){
			$line = substr($line, 0, $len);
			print $fh $line."\n";
		}

		else {
			print $fh $line."\n";
		}

		# Stop if maximum number of reads is met		
		last if $counter == $total;

	}

	if ($format eq 'gzip'){ binmode $fh, ":gzip(none)"; }

}