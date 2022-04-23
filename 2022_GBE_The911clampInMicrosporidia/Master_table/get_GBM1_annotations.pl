#!/usr/bin/perl
## Pombert lab, 2022

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage =<<"USAGE";
$0 -m microsporida_db_file -n ncbi_file -o output.tsv
USAGE
die "\n$usage\n" unless @ARGV;

my $ncbi;
my $msdb;
my $outfile;
GetOptions(
	'm|msdb=s' => \$msdb,
	'n|ncbi=s' => \$ncbi,
	'o|outfile=s' => \$outfile
);

if($msdb){
	open GFF, "<", $msdb or die "Can't open $msdb: $!\n";
	open OUT, ">", "GBM1_MSDB_annotations.tsv" or die "Can't create GBM1_MSDB_annotations.tsv: $!\n";

	print OUT "### Locus tag\tProduct description\n";
	while (my $line = <GFF>){

		chomp $line;
		if ($line =~ /^#/){ next; }
		my @columns = split("\t", $line);
		my $feature = $columns[2];

		if ($feature eq 'protein_coding_gene'){

			my @desc = split(";", $columns[-1]);
			my $locus_tag;
			my $product;

			while (my $tag = shift @desc){
				if ($tag =~ /^ID=(.*)$/){
					$locus_tag = $1;
					$locus_tag =~ s/\D$//;
				}
				if ($tag =~ /^description=(.*) \[/){
					$product = $1;
				}
				elsif ($tag =~ /^description=(.*)/){
					$product = $1;
				}
			}

			print OUT "$locus_tag\t$product\n";

		}
	}
}

if($ncbi){
	open GFF, "<", $ncbi or die "Can't open $ncbi: $!\n";
	open OUT, ">", "GBM1_NCBI_annotations.tsv" or die "Can't create GBM1_NCBI_annotations.tsv: $!\n";

	print OUT "### Locus tag\tProduct description\n";
	while (my $line = <GFF>){

		chomp $line;
		if ($line =~ /^#/){ next; }
		my @columns = split("\t", $line);
		my $feature = $columns[2];

		if ($feature eq 'CDS'){

			my @desc = split(";", $columns[-1]);
			my $locus_tag;
			my $product;

			while (my $tag = shift @desc){
				if ($tag =~ /^locus_tag=(.*)$/){
					$locus_tag = $1;
				}
				if ($tag =~ /^product=(.*)$/){
					$product = $1;
				}
			}

			print OUT "$locus_tag\t$product\n";

		}
	}
}