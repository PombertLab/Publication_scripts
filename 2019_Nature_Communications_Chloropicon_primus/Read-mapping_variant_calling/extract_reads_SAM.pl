#!/usr/bin/perl
## Pombert Lab, 2019
my $version = 0.1;
my $name = 'extract_read_SAM.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $options = <<"OPTIONS";
NAME		$name
VERSION		$version
SYNOPSIS	Create FASTQ datasets for reads mapping/not mapping to a reference based on the SAM alignment
USAGE		parse_reads.pl -s sam -i input.fastq -o output

OPTIONS:
-s (--sam)	SAM alignment
-f (--fastq)	Original FASTQ file(s) used to create the SAM alignment 
-m (--min)	Minimum read length to keep [Default: 1]
OPTIONS
die "\n$options\n" unless @ARGV;
my $sam;
my @fastq;
my $min = 1;
GetOptions(
	's|sam=s' => \$sam,
	'f|fastq=s@{1,}' => \@fastq,
	'm|min=i' => \$min
);
## Init read databases
my %mapped; my %unmapped;
## Working on SAM file
open SAM, "<$sam";
while (my $line = <SAM>){
	chomp $line;
	if ($line =~ /^@/){next;} ## SAM header; skip
	else{
		my @columns = split("\t", $line);
		my $qname = $columns[0];	# QNAME String [!-?A-~]{1,254} Query template NAME
		my $flag = $columns[1];	# FLAG Int [0, 216 − 1] bitwise FLAG
		if ($flag == 4){$unmapped{$qname};}
		else{$mapped{$qname};}
	}
}
## Working on FASTQ file(s)
while (my $fastq = shift@fastq){
	open FASTQ, "<$fastq";
	$fastq =~ s/.fastq$//; $fastq =~ s/.fq$//;
	open MAP, ">$fastq.mapped.fastq";
	open UNMAP, ">$fastq.unmapped.fastq";
	my $count = 0;
	my $read; my $seq; my $qual;
	while (my $line = <FASTQ>){
		chomp $line;
		if ($count == 0){
			$count++;
			if ($line =~ /^@(\S+)/){$read = $1};
		}
		elsif($count == 1){
			$count++;
			$seq = $line;
		} 
		elsif($count == 2){$count++;}
		elsif($count == 3){
			$qual = $line;
			if (length$seq >= $min){
				if (exists $mapped{$read}){
					print MAP '@'."$read\n";
					print MAP "$seq\n";
					print MAP "+\n";
					print MAP "$qual\n";
				}
				elsif (exists $unmapped{$read}){
					print UNMAP '@'."$read\n";
					print UNMAP "$seq\n";
					print UNMAP "+\n";
					print UNMAP "$qual\n";
				}
				else {print "READ $read not present in SAM file\n";}
			}
			$count = 0;
		}
	}
}