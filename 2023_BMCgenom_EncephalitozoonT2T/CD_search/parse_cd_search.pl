#!/usr/bin/perl
my $name = 'parse_cd_search.pl';
my $version = '0.1';
my $updated = '2022-10-25';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = << "USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Parses the content of batch CD-searches by regular expression 
		matching specified terms
		- https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi

REQS		ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz

COMMAND		${name} \
		  -t cddid.tbl \
		  -h *_hitdata.txt \
		  -r methyl mbd \
		  -e 1e-05

OPTIONS:
-t (--tbl)		cddid.tbl file from CD-search with db entries and their description
-h (--hits)		Tab-delimited matches generated by CD-search
-r (--reg)		Regular expression(s) to search for
-e (--eval)		E-value cutoff [Default: 1e-01]
USAGE
die "\n$usage\n" unless @ARGV;

my $cddid;
my @hits;
my @regs;
my $ev_cutoff = 1e-01;
GetOptions(
	't|tbl=s' => \$cddid,
	'h|hits=s@{1,}' => \@hits,
	'r|reg=s@{1,}' => \@regs,
	'e|eval=s' => \$ev_cutoff
);

########## Creating database of cdd IDs and their descriptions #################

open CDD, "<", $cddid or die "Can't open $cddid: $!\n";

my %cdd_ids;
while (my $line = <CDD>){

	chomp $line;
	my @cols = split("\t", $line);
	# 214330	CHL00001	rpoB	RNA polymerase beta subunit	1070

	my $id = $cols[0];
	my $short_name = $cols[1];
	my $desc = $cols[3];
	$cdd_ids{$id} = $desc;

}

########## Parsing CD-search output files ######################################

while (my $hit_file = shift@hits){

	open HIT, "<", $hit_file or die "Can't open $hit_file: $!\n";
	my ($prefix) = $hit_file =~ /^(.*)\.\w+$/;
	my $outfile = $prefix.'.parsed.tsv';
	open OUT, ">", $outfile or die "Can't create $outfile: $!\n";

	while (my $line = <HIT>){
		chomp $line;
		if ($line =~ /^(#|Query)/){
			next;
		}
		elsif ($line =~ /^Q\#/){
			# Query	Hit type	PSSM-ID	From	To	E-Value	Bitscore	Accession	Short name	Incomplete	Superfamily
			my @data = split("\t", $line);
			my $query = $data[0];
			my $hit_type = $data[1];
			my $PSSM_ID= $data[2];
			my $from = $data[3];
			my $to = $data[4];
			my $eval = $data[5];
			my $bitscore = $data[6];
			my $accession = $data[7];
			my $short_name = $data[8];
			my $incomplete = $data[9];
			my $superfamily = $data[10];

			my $description = $cdd_ids{$PSSM_ID};
			if ($eval <= $ev_cutoff){
				foreach my $regex (@regs){
					if ( ($description =~ /$regex/i) or ($short_name =~ /$regex/i) ){
						print OUT $line."\t".$description."\n";
						last;
					}
				}
			}

		}
	}

}