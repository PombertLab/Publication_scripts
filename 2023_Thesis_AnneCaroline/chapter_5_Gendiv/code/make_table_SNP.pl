#!/usr/bin/env perl
## Pombert Lab, Illinois Tech, 2023

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;

my $name = 'make_table_SNP.pl';
my $version = '0.3d';
my $updated = '2023-05-01';

my $usage = << "OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}

COMMAND		${name} \\
		  -l Table_list_SNP.txt \\
		  -acc accessions.tsv \\
		  -ani fastANI.txt \\
		  -mash ./MASH/Mash.mash \\
		  -mummer Mummer/*.report \\
		  -s150 SSRG/SSRG_150/minimap2.varscan2.stats/*.stats \\
		  -s250 SSRG/SSRG_250/minimap2.varscan2.stats/*.stats \\
		  -sill SSRG/SSRG_ILLUMINA/minimap2.varscan2.stats/*.stats \\
		  -e150 snpEFF/VCFs_150/*.csv \\
		  -e250 snpEFF/VCFs_250/*.csv \\
		  -eill snpEFF/VCFs_ILLUMINA/*.csv \\
		  -out table_SNP.tsv
OPTIONS
die "\n$usage\n" unless @ARGV;

my $table_order = 'Table_list_SNP.txt';
my $accessions = 'accessions.tsv';
my $fastani = 'fastANI.txt';
my $mash = './MASH/Mash.mash';
my @mummer_files;
my @stat_files_150;
my @stat_files_250;
my @stat_files_illumina;
my @snpeff_files_150;
my @snpeff_files_250;
my @snpeff_files_illumina;
my $tsv = 'table_SNP.tsv';

GetOptions(
	'l|list=s' => \$table_order,
	'acc=s' => \$accessions,
	'ani=s' => \$fastani,
	'mash=s' => \$mash,
	'm|mummer=s@{1,}' => \@mummer_files,
	's150=s@{1,}' => \@stat_files_150,
	's250=s@{1,}' => \@stat_files_250,
	'sill=s@{1,}' => \@stat_files_illumina,
	'e150=s@{1,}' => \@snpeff_files_150,
	'e250=s@{1,}' => \@snpeff_files_250,
	'eill=s@{1,}' => \@snpeff_files_illumina,
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
		$acc{$isolate}{'gc'} = $gc;
		$acc{$isolate}{'genome'} = $genome;
		$acc{$isolate}{'sra'} = $sra;

	}

}

#####  Mash ##### 
open MASH, "<", $mash or die $!;
while (my $line = <MASH>){

	chomp $line;
	$line =~ s/\.fasta//g;

	my @col = split("\t", $line);
	my ($ref,$path1) = fileparse($col[0]);
	my ($query,$path2) = fileparse($col[1]);
	my $dist = $col[2];
	my $kmer = $col[-1];

	$data{$ref}{$query}{'mash_dist'} = $dist;
	$data{$ref}{$query}{'mash_kmer'} = $kmer;

}

#####  fastANI ##### 
open ANI, "<", $fastani or die $!;
while (my $line = <ANI>){

	chomp $line;
	$line =~ s/\.fasta//g;

	my @col = split("\t", $line);
	my ($ref,$path1) = fileparse($col[0]);
	my ($query,$path2) = fileparse($col[1]);
	my $ani = $col[2];

	$data{$ref}{$query}{'ani'} = $ani;

}

#####  Mummer ##### 
while (my $mummer = shift@mummer_files){
	
	open MUM, "<", $mummer or die $!;

	my ($basename,$path) = fileparse($mummer);
	my ($query,$ref) = $basename =~ /^(\w+)\.vs\.(\w+)\.report/; 

	while (my $line = <MUM>){

		if ($line =~ /^AlignedBases\s+(\d+)\((\S+)\)/){
			$data{$ref}{$query}{'mum_align'} = comma($1);
			$data{$ref}{$query}{'mum_align_pc'} = comma($2);
		}
		elsif ($line =~ /^UnalignedBases\s+(\d+)\((\S+)\)/){
			$data{$ref}{$query}{'mum_unalign'} = comma($1);
			$data{$ref}{$query}{'mum_unalign_pc'} = comma($2);
		}
		elsif ($line =~ /^AvgIdentity\s+(\S+)/){
			$data{$ref}{$query}{'mum_id'} = $1;
		}
		elsif ($line =~ /^TotalSNPs\s+(\d+)/){
			$data{$ref}{$query}{'mum_snps'} = comma($1);
		}
		elsif ($line =~ /^Insertions\s+(\d+)/){
			$data{$ref}{$query}{'mum_insertions'} = comma($1);
		}
		elsif ($line =~ /^Breakpoints\s+(\d+)/){
			$data{$ref}{$query}{'mum_breakpoints'} = comma($1);
		}
		elsif ($line =~ /^Relocations\s+(\d+)/){
			$data{$ref}{$query}{'mum_relocations'} = comma($1);
		}
		elsif ($line =~ /^Translocations\s+(\d+)/){
			$data{$ref}{$query}{'mum_translocations'} = comma($1);
		}
		elsif ($line =~ /^Inversions\s+(\d+)/){
			$data{$ref}{$query}{'mum_inversions'} = comma($1);
		}

	}


}

#####  get_SNPs.pl stats ##### 
while (my $stat = shift@stat_files_150){
	
	open STAT, "<", $stat or die $!;

	my ($basename,$path) = fileparse($stat);
	my ($query,$ref) = $basename =~ /^(\w+).150.R1.fastq.gz.(\w+).*stats/;

	while (my $line = <STAT>){

		chomp $line;

		if ($line =~ /Number of bases covered by at least one read\t(\d+)/){
			$data{$ref}{$query}{'rm_align_150'} = comma($1);
		}
		elsif ($line =~ /Number of bases without coverage\t(\d+)/){
			$data{$ref}{$query}{'rm_unalign_150'} = comma($1);
		}
		elsif ($line =~ /Sequencing breadth \(percentage of bases covered by at least one read\)\t(\S+)/){
			$data{$ref}{$query}{'rm_align_pc_150'} = comma($1);
		}
		elsif ($line =~ /Total number of SNPs \+ indels found:\s+(\d+)/){
			$data{$ref}{$query}{'rm_snps_150'} = comma($1);
		}
		elsif ($line =~ /Average number of SNPs \+ indels per Kb:\s+(\S+)/){
			$data{$ref}{$query}{'rm_snp_perkb_150'} = $1;
		}

	}

}

while (my $stat = shift@stat_files_250){
	
	open STAT, "<", $stat or die $!;

	my ($basename,$path) = fileparse($stat);
	my ($query,$ref) = $basename =~ /^(\w+).250.R1.fastq.gz.(\w+).*stats/;

	while (my $line = <STAT>){

		chomp $line;

		if ($line =~ /Number of bases covered by at least one read\t(\d+)/){
			$data{$ref}{$query}{'rm_align_250'} = comma($1);
		}
		elsif ($line =~ /Number of bases without coverage\t(\d+)/){
			$data{$ref}{$query}{'rm_unalign_250'} = comma($1);
		}
		elsif ($line =~ /Sequencing breadth \(percentage of bases covered by at least one read\)\t(\S+)/){
			$data{$ref}{$query}{'rm_align_pc_250'} = comma($1);
		}
		elsif ($line =~ /Total number of SNPs \+ indels found:\s+(\d+)/){
			$data{$ref}{$query}{'rm_snps_250'} = comma($1);
		}
		elsif ($line =~ /Average number of SNPs \+ indels per Kb:\s+(\S+)/){
			$data{$ref}{$query}{'rm_snp_perkb_250'} = $1;
		}

	}

}

while (my $stat = shift@stat_files_illumina){
	
	open STAT, "<", $stat or die $!;

	my ($basename,$path) = fileparse($stat);
	my ($query,$ref) = $basename =~ /^(\w+)\..*.fastq.gz.(\w+).*stats/;

	while (my $line = <STAT>){

		chomp $line;

		if ($line =~ /Number of bases covered by at least one read\t(\d+)/){
			$data{$ref}{$query}{'rm_align_illumina'} = comma($1);
		}
		elsif ($line =~ /Number of bases without coverage\t(\d+)/){
			$data{$ref}{$query}{'rm_unalign_illumina'} = comma($1);
		}
		elsif ($line =~ /Sequencing breadth \(percentage of bases covered by at least one read\)\t(\S+)/){
			$data{$ref}{$query}{'rm_align_pc_illumina'} = comma($1);
		}
		elsif ($line =~ /Total number of SNPs \+ indels found:\s+(\d+)/){
			$data{$ref}{$query}{'rm_snps_illumina'} = comma($1);
		}
		elsif ($line =~ /Average number of SNPs \+ indels per Kb:\s+(\S+)/){
			$data{$ref}{$query}{'rm_snp_perkb_illumina'} = $1;
		}

	}

}

##### snpEff ##### 

while (my $snpeff = shift@snpeff_files_150){
	
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
			$data{$ref}{$query}{'rm_snp_del_150'} = $del;
		}
		elsif ($line =~ /^INS , (\d+)/){
			$insert = comma($1);
			$data{$ref}{$query}{'rm_snp_ins_150'} = $insert;
		}
		elsif ($line =~ /^MISSENSE , (\d+)/){
			$missense = comma($1);
			$data{$ref}{$query}{'rm_snp_missense_150'} = $missense;
		}
		elsif ($line =~ /^NONSENSE , (\d+)/){
			$nonsense = comma($1);
			$data{$ref}{$query}{'rm_snp_nonsense_150'} = $nonsense;
		}
		elsif ($line =~ /^SILENT , (\d+)/){
			$silent = comma($1);
			$data{$ref}{$query}{'rm_snp_silent_150'} = $silent;
		}


		if (!defined $del){
			$data{$ref}{$query}{'rm_snp_del_150'} = 0;
		}
		if (!defined $insert){
			$data{$ref}{$query}{'rm_snp_ins_150'} = 0;
		}
		if (!defined $missense){
			$data{$ref}{$query}{'rm_snp_missense_150'} = 0;
		}
		if (!defined $nonsense){
			$data{$ref}{$query}{'rm_snp_nonsense_150'} = 0;
		}
		if (!defined $silent){
			$data{$ref}{$query}{'rm_snp_silent_150'} = 0;
		}

		## 
		my $syn = $data{$ref}{$query}{'rm_snp_silent_150'};
		my $nonsyn = $data{$ref}{$query}{'rm_snp_missense_150'};
		my $ratio;

		if ( ($syn and $nonsyn) == 0 ){
			$ratio = 0;
		}
		elsif ($nonsyn == 0){
			$ratio = '---';
		}
		else {
			$ratio = sprintf("%.3f", $syn/$nonsyn);
		}

		$data{$ref}{$query}{'rm_syn_ratio_150'} = $ratio;

	}

}

while (my $snpeff = shift@snpeff_files_250){
	
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
			$data{$ref}{$query}{'rm_snp_del_250'} = $del;
		}
		elsif ($line =~ /^INS , (\d+)/){
			$insert = comma($1);
			$data{$ref}{$query}{'rm_snp_ins_250'} = $insert;
		}
		elsif ($line =~ /^MISSENSE , (\d+)/){
			$missense = comma($1);
			$data{$ref}{$query}{'rm_snp_missense_250'} = $missense;
		}
		elsif ($line =~ /^NONSENSE , (\d+)/){
			$nonsense = comma($1);
			$data{$ref}{$query}{'rm_snp_nonsense_250'} = $nonsense;
		}
		elsif ($line =~ /^SILENT , (\d+)/){
			$silent = comma($1);
			$data{$ref}{$query}{'rm_snp_silent_250'} = $silent;
		}


		if (!defined $del){
			$data{$ref}{$query}{'rm_snp_del_250'} = 0;
		}
		if (!defined $insert){
			$data{$ref}{$query}{'rm_snp_ins_250'} = 0;
		}
		if (!defined $missense){
			$data{$ref}{$query}{'rm_snp_missense_250'} = 0;
		}
		if (!defined $nonsense){
			$data{$ref}{$query}{'rm_snp_nonsense_250'} = 0;
		}
		if (!defined $silent){
			$data{$ref}{$query}{'rm_snp_silent_250'} = 0;
		}

		## 
		my $syn = $data{$ref}{$query}{'rm_snp_silent_250'};
		my $nonsyn = $data{$ref}{$query}{'rm_snp_missense_250'};
		my $ratio;

		if ( ($syn and $nonsyn) == 0 ){
			$ratio = 0;
		}
		elsif ($nonsyn == 0){
			$ratio = '---';
		}
		else {
			$ratio = sprintf("%.3f", $syn/$nonsyn);
		}

		$data{$ref}{$query}{'rm_syn_ratio_250'} = $ratio;

	}

}

while (my $snpeff = shift@snpeff_files_illumina){
	
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
			$data{$ref}{$query}{'rm_snp_del_illumina'} = $del;
		}
		elsif ($line =~ /^INS , (\d+)/){
			$insert = comma($1);
			$data{$ref}{$query}{'rm_snp_ins_illumina'} = $insert;
		}
		elsif ($line =~ /^MISSENSE , (\d+)/){
			$missense = comma($1);
			$data{$ref}{$query}{'rm_snp_missense_illumina'} = $missense;
		}
		elsif ($line =~ /^NONSENSE , (\d+)/){
			$nonsense = comma($1);
			$data{$ref}{$query}{'rm_snp_nonsense_illumina'} = $nonsense;
		}
		elsif ($line =~ /^SILENT , (\d+)/){
			$silent = comma($1);
			$data{$ref}{$query}{'rm_snp_silent_illumina'} = $silent;
		}


		if (!defined $del){
			$data{$ref}{$query}{'rm_snp_del_illumina'} = 0;
		}
		if (!defined $insert){
			$data{$ref}{$query}{'rm_snp_ins_illumina'} = 0;
		}
		if (!defined $missense){
			$data{$ref}{$query}{'rm_snp_missense_illumina'} = 0;
		}
		if (!defined $nonsense){
			$data{$ref}{$query}{'rm_snp_nonsense_illumina'} = 0;
		}
		if (!defined $silent){
			$data{$ref}{$query}{'rm_snp_silent_illumina'} = 0;
		}

		## 
		my $syn = $data{$ref}{$query}{'rm_snp_silent_illumina'};
		my $nonsyn = $data{$ref}{$query}{'rm_snp_missense_illumina'};
		my $ratio;

		if ( ($syn and $nonsyn) == 0 ){
			$ratio = 0;
		}
		elsif ($nonsyn == 0){
			$ratio = '---';
		}
		else {
			$ratio = sprintf("%.3f", $syn/$nonsyn);
		}

		$data{$ref}{$query}{'rm_syn_ratio_illumina'} = $ratio;

	}

}

#####  Working on TSV table ##### 

open TSV, ">", $tsv or die $!;

my @header = (
	## From Table_list.txt
	'QUERY',
	'REF',
	## From accessions.txt
	'REF_DESC',
	'QUERY_GENOME',
	'REF_GENOME',
	'QUERY_SRA',
	'REF_SRA',
	'QUERY_GC',
	'REF_GC',
	## Distance
	'ani',
	'mash_dist',
	'mash_kmer',
	## Mummer
	'mum_align',
	'mum_align_pc',
	'mum_unalign_pc',
	'mum_id',
	'mum_snps',
	'mum_insertions',
	'mum_relocations',
	'mum_translocations',
	'mum_inversions',
	'mum_breakpoints',
	## RM 150 bp
	'rm_align_150',
	'rm_unalign_150',
	'rm_align_pc_150',
	'rm_snps_150',
	'rm_snp_perkb_150',
	'rm_snp_del_150',
	'rm_snp_ins_150',
	'rm_snp_silent_150',
	'rm_snp_missense_150',
	'rm_snp_nonsense_150',
	'rm_syn_ratio_150',
	## RM 250 bp
	'rm_align_250',
	'rm_unalign_250',
	'rm_align_pc_250',
	'rm_snps_250',
	'rm_snp_perkb_250',
	'rm_snp_del_250',
	'rm_snp_ins_250',
	'rm_snp_silent_250',
	'rm_snp_missense_250',
	'rm_snp_nonsense_250',
	'rm_syn_ratio_250',
	## RM illumina
	'rm_align_illumina',
	'rm_unalign_illumina',
	'rm_align_pc_illumina',
	'rm_snps_illumina',
	'rm_snp_perkb_illumina',
	'rm_snp_del_illumina',
	'rm_snp_ins_illumina',
	'rm_snp_silent_illumina',
	'rm_snp_missense_illumina',
	'rm_snp_nonsense_illumina',
	'rm_syn_ratio_illumina',
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
		my @col = split("\t", $line);
		my $query = $col[0];
		my $ref = $col[1];

		my $description = ${acc}{$ref}{'description'};
		my $q_genome = ${acc}{$query}{'genome'};
		my $r_genome = ${acc}{$ref}{'genome'};
		my $q_sra = ${acc}{$query}{'sra'};
		my $r_sra = ${acc}{$ref}{'sra'};
		my $q_gc = ${acc}{$query}{'gc'};
		my $r_gc = ${acc}{$ref}{'gc'};

		print TSV "$query\t";
		print TSV "$ref\t";
		print TSV "$description\t";
		print TSV "$q_genome\t";
		print TSV "$r_genome\t";
		print TSV "$q_sra\t";
		print TSV "$r_sra\t";
		print TSV "$q_gc\t";
		print TSV "$r_gc\t";

		for (9..$#header-1){
			my $key = $header[$_];
			print TSV "$data{$ref}{$query}{$key}\t";
		}
		print TSV "$data{$ref}{$query}{$header[-1]}\n";
		
	}

}

################ Subroutine(s
# Adds commas to numbers; from the Perl Cookbook (O'Reilly)
sub comma {
	my $string = reverse $_[0];
	$string =~ s/(\d{3})(?=\d)(?!\d*\.)/$1,/g;
	return scalar (reverse ($string));
}