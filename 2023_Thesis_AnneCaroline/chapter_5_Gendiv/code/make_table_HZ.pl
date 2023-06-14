#!/usr/bin/env perl
## Pombert Lab, Illinois Tech, 2023

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;

my $name = 'make_table_SNP.pl';
my $version = '0.1';
my $updated = '2023-05-23';

my $usage = << "OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}

COMMAND		${name} \\
		  -l Table_list_HZ.txt \\
		  -acc accessions.tsv \\
		  -stats RM/HZ/minimap2.varscan2.stats/*.stats \\
		  -eff snpEFF/VCFs_HZ/*.csv \\
		  -out table_HZ.tsv
OPTIONS
die "\n$usage\n" unless @ARGV;

my $table_order = 'Table_list_HZ.txt';
my $accessions = 'accessions.tsv';
my @stat_files;
my @snpeff_files;
my $tsv = 'table_HZ.tsv';

GetOptions(
	'l|list=s' => \$table_order,
	'acc=s' => \$accessions,
	'stats=s@{1,}' => \@stat_files,
	'eff=s@{1,}' => \@snpeff_files,
	'o|out=s' => \$tsv
);

##### databases #####
my %acc;
my %data;

##### Accessions #####
open ACC, "<", $accessions or die $!;
while (my $line = <ACC>){

	chomp $line;
	if ($line !~ /^#/){

		my @col = split("\t", $line);
		my $isolate = $col[0];
		my $genotype = $col[1];
		my $description = $col[2];
		my $gc = $col[3];
		my $genome = $col[4];
		my $sra = $col[5];

		$acc{$isolate}{'description'} = $description;
		$acc{$isolate}{'genotype'} = $genotype;
		$acc{$isolate}{'genome'} = $genome;
		$acc{$isolate}{'sra'} = $sra;
		$acc{$isolate}{'gc'} = $gc;

	}

}

#####  get_SNPs.pl stats ##### 
while (my $stat = shift@stat_files){
	
	open STAT, "<", $stat or die $!;

	my ($basename,$path) = fileparse($stat);
	my ($query,$ref) = $basename =~ /^(\w+)\..*.fastq.gz.(\w+).*stats/;

	while (my $line = <STAT>){

		chomp $line;

		if ($line =~ /Number of bases covered by at least one read\t(\d+)/){
			$data{$ref}{$query}{'rm_align'} = comma($1);
		}
		elsif ($line =~ /Number of bases without coverage\t(\d+)/){
			$data{$ref}{$query}{'rm_unalign'} = comma($1);
		}
		elsif ($line =~ /Sequencing breadth \(percentage of bases covered by at least one read\)\t(\S+)/){
			$data{$ref}{$query}{'rm_align_pc'} = comma($1);
		}
		elsif ($line =~ /Total number of SNPs \+ indels found:\s+(\d+)/){
			$data{$ref}{$query}{'rm_snps'} = comma($1);
		}
		elsif ($line =~ /Average number of SNPs \+ indels per Kb:\s+(\S+)/){
			$data{$ref}{$query}{'rm_snp_perkb'} = $1;
		}

	}

}

##### snpEff ##### 

while (my $snpeff = shift@snpeff_files){
	
	open EFF, "<", $snpeff or die $!;

	my ($basename,$path) = fileparse($snpeff);
	my ($query,$ref) = $basename =~ /^(\w+)\.vs\.(\w+)\.summary\.csv/; 

	my $del = undef;
	my $insert = undef;
	my $missense = undef;
	my $nonsense = undef;
	my $silent = undef;

	while (my $line = <EFF>){

		chomp $line;

		if ($line =~ /^DEL , (\d+)/){
			$del = comma($1); 
			$data{$ref}{$query}{'rm_snp_del'} = $del;
		}
		elsif ($line =~ /^INS , (\d+)/){
			$insert = comma($1);
			$data{$ref}{$query}{'rm_snp_ins'} = $insert;
		}
		elsif ($line =~ /^MISSENSE , (\d+)/){
			$missense = comma($1);
			$data{$ref}{$query}{'rm_snp_missense'} = $missense;
		}
		elsif ($line =~ /^NONSENSE , (\d+)/){
			$nonsense = comma($1);
			$data{$ref}{$query}{'rm_snp_nonsense'} = $nonsense;
		}
		elsif ($line =~ /^SILENT , (\d+)/){
			$silent = comma($1);
			$data{$ref}{$query}{'rm_snp_silent'} = $silent;
		}


		if (!defined $del){
			$data{$ref}{$query}{'rm_snp_del'} = 0;
		}
		if (!defined $insert){
			$data{$ref}{$query}{'rm_snp_ins'} = 0;
		}
		if (!defined $missense){
			$data{$ref}{$query}{'rm_snp_missense'} = 0;
		}
		if (!defined $nonsense){
			$data{$ref}{$query}{'rm_snp_nonsense'} = 0;
		}
		if (!defined $silent){
			$data{$ref}{$query}{'rm_snp_silent'} = 0;
		}

	}

}

#####  Working on TSV table ##### 

open TSV, ">", $tsv or die $!;

my @header = (
	## From Table_list_HZ.txt and accessions.tsv
	'ISOLATE',
	'ISO_DESC',
	'ISO_ACCESSION',
	'ISO_SRA',
	## RM illumina
	'rm_align',
	'rm_unalign',
	'rm_align_pc',
	'rm_snps',
	'rm_snp_perkb',
	'rm_snp_del',
	'rm_snp_ins',
	'rm_snp_silent',
	'rm_snp_missense',
	'rm_snp_nonsense',
);

for (0..$#header-1){
	print TSV "$header[$_]\t";
}
print TSV "$header[-1]\n";

open LIST, "<", $table_order or die $!;
while (my $line = <LIST>){
	
	chomp $line;
	unless ($line =~ /^#/){

		chomp $line;
		my $isolate = $line;
		my $iso_desc = $acc{$isolate}{'description'};
		my $iso_accession = $acc{$isolate}{'genome'};
		my $iso_sra = $acc{$isolate}{'sra'};

		print TSV "$isolate\t";
		print TSV "$iso_desc\t";
		print TSV "$iso_accession\t";
		print TSV "$iso_sra\t";
		for (4..$#header-1){
			my $key = $header[$_];
			print TSV "$data{$isolate}{$isolate}{$key}\t";
		}
		print TSV "$data{$isolate}{$isolate}{$header[-1]}\n";
		
	}

}

################ Subroutine(s
# Adds commas to numbers; from the Perl Cookbook (O'Reilly)
sub comma {
	my $string = reverse $_[0];
	$string =~ s/(\d{3})(?=\d)(?!\d*\.)/$1,/g;
	return scalar (reverse ($string));
}