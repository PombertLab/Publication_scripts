#!/usr/bin/env perl
## Pombert Lab, Illinois Tech, 2023
my $name = 'snpeff_metrics.pl';
my $version = '0.1';
my $updated = '2023-04-27';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Parses the output of SnpEff annotated VCF files to calculate
		transitions, transversions and other metrics

COMMAND		${name} \\
		  -v *.annotations.vcf \\
		  -o outdir \\
		  -p prefix \\
		  -l Table_list_SNP.txt \\
		  -s locus_tag_list.txt \\
		  -v

-v (--vcf)	Annotated snpEFF VCF file(s) to parse 
-o (--out)	Output directory [Default: ./]
-p (--prefix)	Output files prefix [Default: metrics]
-l (--list)	Table output order [Default: Table_list_SNP.txt]
-s (--strand)	Check for strandedness from tab-delimited list: locus_tag and +/-
--verbose	Turn on debugging statements
OPTIONS
die "\n$usage\n" unless @ARGV;

my @vcf_files;
my $outdir = './';
my $prefix = 'metrics';
my $table_order = 'Table_list_SNP.txt';
my $strands;
my $verbose;
GetOptions(
	'v|vcf=s@{1,}' => \@vcf_files,
	'o|out=s' => \$outdir,
	'p|prefix=s' => \$prefix,
	'l|list=s' => \$table_order,
	's|strand=s' => \$strands,
	'verbose' => \$verbose
);

############## strandeness
my %loci;

if ($strands){

	open LOCI, '<', $strands or die $!;
	while (my $line = <LOCI>){

		chomp $line;
		my ($locus,$strand) = split("\t", $line); 

		$loci{$locus} = $strand;

	}

}

############## Outputs

unless (-d $outdir){
	mkdir ($outdir, 0755) or die $!;
}

my $table = $outdir.'/'.$prefix.'.table.tsv';
my $data = $outdir.'/'.$prefix.'.data.txt';
my $data_adj = $outdir.'/'.$prefix.'.data_adj.txt';

open TSV, '>', $table or die $!;
open DT, '>', $data or die $!;
open ADJ, '>', $data_adj or die $!;

print TSV 'QUERY'."\t".'REFERENCE';
print TSV "\t".'SNP_SUM';
print TSV "\t".'Ts'."\t".'Ts %';
print TSV "\t".'Tv'."\t".'Tv %';
print TSV "\t".'Ts/Tv';

print TSV "\t".'A > C'."\t";
print TSV "\t".'A > G'."\t";
print TSV "\t".'A > T'."\t";

print TSV "\t".'C > A'."\t";
print TSV "\t".'C > G'."\t";
print TSV "\t".'C > T'."\t";

print TSV "\t".'G > A'."\t";
print TSV "\t".'G > C'."\t";
print TSV "\t".'G > T'."\t";

print TSV "\t".'T > A'."\t";
print TSV "\t".'T > C'."\t";
print TSV "\t".'T > G'."\n";

############## Working on VCF files

my %db;			## Regular db

my %db_strand;	## Stores data based on ref orientation
				## to check for possible strand skews
				## Mutation could occur on either strand
				## but would an effect be noticeable based
				## on strandeness?
my %dnds;

while (my $vcf_file = shift @vcf_files){
	
	open VCF, '<', $vcf_file or die $!;
	my ($vcf, $path) = fileparse($vcf_file);

	## Making sure that the keys are generated even if VCFs are empty
	## (if genomes are identical)
	foreach my $base ('A','T','G','C'){
		foreach ('A','T','G','C'){
			$db{$vcf}{$base}{$_} = 0;
		}
	}

	foreach my $base ('A','T','G','C'){
		foreach ('A','T','G','C'){
			$db_strand{$vcf}{$base}{$_} = 0;
		}
	}

	## Working on VCF content
	while (my $line = <VCF>){
		
		unless ($line =~ /^#/){
			
			my @cols = split("\t", $line);
			my $locus  = $cols[0];
			my $position = $cols[1];
			my $id = $cols[2];
			my $ref = $cols[3];
			my $alt = $cols[4];
			my $qual = $cols[5];
			my $filer = $cols[6];
			my $info = $cols[7];

			## Checking for INDEL
			if ( (length($ref) > 1) or (length($alt) > 1) ){
				next;
			}

			## Checking for single nucleotide variants
			else {

				$db{$vcf}{$ref}{$alt} += 1;
				
				if ($info =~ /missense_variant\|\w+\|(\w+)/){

					my $locus = $1;
					$dnds{$vcf}{'dn'} += 1;

					if ($strands){
						if (exists $loci{$locus}){
							if ($loci{$locus} eq '+'){
								$db_strand{$vcf}{$ref}{$alt} += 1;
							}
							else {
								$ref =~ tr/ATGCatgc/TACGtacg/;
								$alt =~ tr/ATGCatgc/TACGtacg/;
								$db_strand{$vcf}{$ref}{$alt} += 1;
							}
						}
						else {
							if ($verbose){
								print "Locus $locus not found\n";
							}
						}
					}


				}
				elsif ($info =~ /synonymous_variant\|\w+\|(\w+)/){
					
					my $locus = $1;
					$dnds{$vcf}{'ds'} += 1;	

					if ($strands){
						if (exists $loci{$locus}){
							if ($loci{$locus} eq '+'){
								$db_strand{$vcf}{$ref}{$alt} += 1;
							}
							else {
								$ref =~ tr/ATGCatgc/TACGtacg/;
								$alt =~ tr/ATGCatgc/TACGtacg/;
								$db_strand{$vcf}{$ref}{$alt} += 1;
							}
						}
						else {
							if ($verbose){
								print "Locus $locus not found\n";
							}
						}
					}

				}
				elsif ($info =~ /stop_gained/){
					$dnds{$vcf}{'stop'} += 1;	
				}
				elsif ($info =~ /intergenic/){
					$dnds{$vcf}{'intergenic'} += 1;	
				}
				elsif ($info =~ /(up|down)stream_gene_variant/){
					$dnds{$vcf}{'updown_gene_variant'} += 1;	
				}

				## Transitions and transversions				
				if ( ($ref =~ /[AGag]/) and ($alt =~ /[AGag]/) ){
					$dnds{$vcf}{'transitions'} += 1;
				}
				elsif ( ($ref =~ /[CTct]/) and ($alt =~ /[CTct]/) ){
					$dnds{$vcf}{'transitions'} += 1;
				}
				else {
					$dnds{$vcf}{'transversions'} += 1;
				}

			}

		}

	}

}

## Print summary
foreach my $vcf (sort (keys %db)){

	## Check for upstream or downtream variants
	unless ($dnds{$vcf}{'updown_gene_variant'}){
		$dnds{$vcf}{'updown_gene_variant'} = 0;
	}
	unless ($dnds{$vcf}{'stop'}){
		$dnds{$vcf}{'stop'} = 0;
	}

	## Calculating dn/ds ratios (avoiding div by zero)
	my $dnds;
	if ($dnds{$vcf}{'ds'}){
		$dnds = sprintf("%.3f", $dnds{$vcf}{'dn'} / $dnds{$vcf}{'ds'});
	}
	else{
		$dnds ='N/A (div0)';
		$dnds{$vcf}{'ds'} = 0;
	}

	unless ($dnds{$vcf}{'dn'}){
		$dnds{$vcf}{'dn'} = 0;
	}

	## Calculating total number of SNPs
	my $snp_sum = 0;
	foreach my $ref_base (keys %{$db{$vcf}}){
		foreach my $alt_base (keys %{$db{$vcf}{$ref_base}}){
			$snp_sum += $db{$vcf}{$ref_base}{$alt_base};
		}
	}

	my $snp_sum_adj = 0;
	if ($strands){
		foreach my $ref_base (keys %{$db_strand{$vcf}}){
			foreach my $alt_base (keys %{$db_strand{$vcf}{$ref_base}}){
				$snp_sum_adj += $db_strand{$vcf}{$ref_base}{$alt_base};
			}
		}
	}

	# Print to summary

	print DT '####################   '.$vcf.'   ####################'."\n\n";
	print DT $vcf."\t".'SNP total: '.comma($snp_sum)."\n";

	unless ($dnds{$vcf}{'transitions'}){
		$dnds{$vcf}{'transitions'} = 0;
	}
	unless ($dnds{$vcf}{'transversions'}){
		$dnds{$vcf}{'transversions'} = 0;
	}

	my $tstv = 0;
	unless ($dnds{$vcf}{'transversions'} == 0){
		$tstv = sprintf("%.3f", $dnds{$vcf}{'transitions'} / $dnds{$vcf}{'transversions'});
	}

	print DT $vcf."\t".'Transitions (Ts): '.comma($dnds{$vcf}{'transitions'})."\n";
	print DT $vcf."\t".'Transversions (Tv): '.comma($dnds{$vcf}{'transversions'})."\n";
	print DT $vcf."\t".'Ts/Tv ratio: '.comma($tstv)."\n";
	print DT $vcf."\t".'Intergenic: '.comma($dnds{$vcf}{'intergenic'})."\n";
	print DT $vcf."\t".'Upstream/downstream gene variant: '.comma($dnds{$vcf}{'updown_gene_variant'})."\n";
	print DT $vcf."\t".'Synonymous: '.comma($dnds{$vcf}{'ds'})."\n";
	print DT $vcf."\t".'Non-synonymous: '.comma($dnds{$vcf}{'dn'})."\n";
	print DT $vcf."\t".'Stop codons: '.comma($dnds{$vcf}{'stop'})."\n";
	print DT $vcf."\t".'Global dN/dS: '.comma($dnds)."\n\n";

	#### Calculating mutation biases
	# Summary, might include ambiguous bases
	foreach my $ref_base (sort (keys %{$db{$vcf}}) ){

		foreach my $alt_base (sort (keys %{$db{$vcf}{$ref_base}}) ){

			my $mut_count = 0;
			if ($db{$vcf}{$ref_base}{$alt_base}){
				$mut_count = $db{$vcf}{$ref_base}{$alt_base};
			}

			my $mut_ratio = 0;
			if ($snp_sum != 0){
				$mut_ratio = ($mut_count/$snp_sum)*100;
			}
			$mut_ratio = sprintf("%.2f %%", $mut_ratio);

			$mut_count = comma($mut_count);
			print DT $vcf."\t".$ref_base.' > '.$alt_base."\t".$mut_count."\t".$mut_ratio."\n";

		}

	}

	print DT "\n";

	## Checking metrics adjusted for strandedness
	if ($strands){
	
		foreach my $ref_base (sort (keys %{$db_strand{$vcf}}) ){

			foreach my $alt_base (sort (keys %{$db_strand{$vcf}{$ref_base}}) ){

				my $mut_count = 0;
				if ($db_strand{$vcf}{$ref_base}{$alt_base}){
					$mut_count = $db_strand{$vcf}{$ref_base}{$alt_base};
				}

				my $mut_ratio = 0;
				if ($snp_sum_adj != 0){
					$mut_ratio = ($mut_count/$snp_sum_adj)*100;
				}
				$mut_ratio = sprintf("%.2f %%", $mut_ratio);

				$mut_count = comma($mut_count);
				print ADJ $vcf."\t".$ref_base.' > '.$alt_base."\t".$mut_count."\t".$mut_ratio."\t".' (strand-adjusted)'."\n";

			}

		}
	
		print ADJ "\n";
	
	}


}

## Print table
open LIST, '<', $table_order or die $!;

while (my $line = <LIST>){
	
	chomp $line;

	next if ($line =~ /^#/); 

	my ($query,$ref) = split("\t",$line);
	my $vcf = "$query.vs.$ref.annotations.vcf";

	print TSV $query."\t".$ref;

	## Calculating total number of SNPs
	my $snp_sum = 0;
	foreach my $ref_base (keys %{$db{$vcf}}){
		foreach my $alt_base (keys %{$db{$vcf}{$ref_base}}){
			$snp_sum += $db{$vcf}{$ref_base}{$alt_base};
		}
	}

	unless ($dnds{$vcf}{'transitions'}){
		$dnds{$vcf}{'transitions'} = 0;
	}
	unless ($dnds{$vcf}{'transversions'}){
		$dnds{$vcf}{'transversions'} = 0;
	}


	my $ts = $dnds{$vcf}{'transitions'};
	my $tv = $dnds{$vcf}{'transversions'};

	my $tstv = 0;
	unless ($tv == 0){
		$tstv = sprintf("%.3f", ($ts/$tv));
	}

	my $ts_ratio = 0;
	my $tv_ratio = 0;
	
	if ($snp_sum == 0){
		$ts_ratio = '---';
		$ts_ratio = '---';
		$tstv = '---';
	}
	else {
		$ts_ratio = sprintf("%.2f %%", ($ts/$snp_sum)*100);
		$tv_ratio = sprintf("%.2f %%", ($tv/$snp_sum)*100);
	}

	print TSV "\t".$snp_sum;
	print TSV "\t".$ts."\t".$ts_ratio;
	print TSV "\t".$tv."\t".$tv_ratio;
	print TSV "\t".$tstv;

	## Table, only ATGC 
	foreach my $ref_base ( sort ('A','T','G','C') ){

		foreach my $alt_base (sort ('A','T','G','C') ){

			if ($ref_base eq $alt_base){
				next;
			}

			else {

				my $mut_count = 0;
				if ($db{$vcf}{$ref_base}{$alt_base}){
					$mut_count = $db{$vcf}{$ref_base}{$alt_base};
				}

				my $mut_ratio = '---';
				if ($snp_sum != 0){
					$mut_ratio = sprintf("%.2f %%", ($mut_count/$snp_sum)*100);
				}

				print TSV "\t".$mut_count."\t".$mut_ratio;

			}

		}

	}

	print TSV "\n";

}

#### Subs
sub comma {
	my $comma = reverse $_[0];
	$comma =~ s/(\d{3})(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $comma;
}