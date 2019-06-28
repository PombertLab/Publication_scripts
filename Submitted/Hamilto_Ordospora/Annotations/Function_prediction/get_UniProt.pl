#!/usr/bin/perl
## Pombert Lab, 2018
my $name = 'get_UniProt.pl';
my $version = 0.1;

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Downloads the SwissProt and/or trEMBL databases from UniProt
EXAMPLE		get_UniProt.pl -s -t -n 20 -l download.log 

OPTIONS:
-s (--swiss)	Download Swiss-Prot
-t (--trembl)	Download trEMBL
-n (--nice)	Linux Process Priority [Default: 20] ## Runs downloads in the background
-l (--log)	Print download information to log file

OPTIONS
die $usage unless @ARGV;

my $nice = 20;
my $swiss;
my $trembl;
my $log;
GetOptions(
	'n|nice=i' => \$nice,
	's|swiss' => \$swiss,
	't|trembl' => \$trembl,
	'l|log=s' => \$log
);

## Creating logs
if ($log){open LOG, ">$log";}

## Downloading SwissProt
if ($swiss){
	print "\nDownloading SwissProt...\n\n";
	if ($log){my $date = `date`; print LOG "Downloading SwissProt on $date";}
	system "nice -n $nice wget -c ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz";
	if ($log){my $size = `du -h uniprot_sprot.fasta.gz`; print LOG "$size";}
}

## Downloading trEMBL
if ($trembl){
	print "\nDownloading trEMBL database. This will take a while...\n\n";
	if ($log){my $date = `date`; print LOG "Downloading trEMBL on $date";	}
	system "nice -n $nice wget -c ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz";
	if ($log){my $size = `du -h uniprot_trembl.fasta.gz`; print LOG "$size";}
}

## Decompressing the downloaded databases
if ($swiss){
	print "\nDecompressing SwissProt db...\n\n";
	system "gunzip uniprot_sprot.fasta.gz";
}
if ($trembl){
	print "\nDecompressing trEMBL db...\n\n";
	system "uniprot_trembl.fasta.gz";
}
## Updating logs
if ($log){
	my $end = `date`;
	print LOG "Finished downloading database(s) on $end";
}
