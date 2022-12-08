#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2022
my $name = 'get_ppiscores.pl';
my $version = '0.3a';
my $updated = '2022-09-29';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Find;
use File::Basename;
use File::Path qw(make_path);
use Cwd qw(abs_path);

my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Get PPIscores from Megadock results

EXAMPLE		${name} \\
		  -m MDOCK \\
		  -o output_dir \\
		  -n 10000 \\
		  -c 10 \\
		  -p Eint_products.tsv \\
		  -t ppiscores.tsv

GENERAL OPTIONS:
-m (--mdock)	MDOCK output directory # From run_megadock.pl
-n (--npred)	Number of output predictions ?? [Default: 10000]
-c (--cutoff)	Minimum PPIscore to keep [Default: 5]
-t (--tsv)	Name of the tsv output file [Default: ppiscores.tsv]
-p (--products)	Tab-delimited list of products for the ligands queried
-v (--verbose)	Adds verbosity
USAGE
die "\n$usage\n" unless @ARGV;

my $megadock;
my $outdir = 'PPIs';
my $tsvfile = 'ppiscores.tsv';
my $npred = 10000;
my $ppicutoff = 5;
my $products;
my $verbose;
GetOptions(
	'm|mdock=s' => \$megadock,
	'o|outdir=s' => \$outdir,
	'n|npred=i' => \$npred,
	't|tsv=s' => \$tsvfile,
	'c|cutoff=i' => \$ppicutoff,
	'p|products=s' => \$products,
	'v|verbose' => \$verbose
);

##### Check if output directory / TSV file can be created
$outdir =~ s/\/$//;
$outdir = abs_path($outdir);
unless (-d $outdir) {
	make_path($outdir, {mode => 0755}) or die "Can't create $outdir: $!\n";
}

##### Create a database of ligands and their products
## to help produce a more informative ouput
my %products;
if ($products){
	open PROD, "<", $products or die "Can't open $products: $!\n";
	while (my $line = <PROD>){
		chomp $line;
		unless ($line =~ /^#/){
			my @data = split("\t", $line);
			my $locus = $data[0];
			my $product = $data[1];
			$products{$locus} = $product;
		}
	}
}

##### Create a hash to store PPIscores and rank them
my %rankings;
$tsvfile = $outdir.'/'.$tsvfile;

# Check if ppiscores.tsv exists; if so preload data (save on computation time)
if (-e $tsvfile){
	open PPI, "<", $tsvfile or die "Can't open $tsvfile: $!\n";
	while (my $line = <PPI>){
		chomp $line;
		if ($line =~ /^#/){ next; }
		else {
			my @data = split("\t", $line);
			my $rec = $data[0];
			my $lig = $data[1];
			my $score = $data[2];
			$rankings{$rec}{$lig} = $score;
		}
	}
	close PPI;
}

##### Getting results from run_megadock.pl output directory
my %mdocks;
find (\&mdock, $megadock);
sub mdock { 
	my $file = $File::Find::name;
	unless (-d $file){
		if ($file =~ /\.out$/){
			my ($receptor) = $file =~ /([A-Za-z0-9_-]+)\.vs/;
			push (@{$mdocks{$receptor}}, $file);
		}
	}
}

# Creating blank ppiscores.tsv 
open TSV, ">", $tsvfile or die "Can't create $tsvfile: $!\n";
print TSV '### Receptor'."\t".'Ligand'."\t".'PPIscore';
if ($products){
	print TSV "\t".'Product';
}
print TSV "\n";

##### Working on megadock files
my @receptors = sort (keys %mdocks);

print "\n### Getting PPIscores\n";

while (my $receptor = shift @receptors){

	if ($verbose){
		print "\n  Working on receptor: $receptor ...\n";
	}

	## Create a hash to store PPIscores and rank them
	my $ligand_count = 0;

	## Iterating through the megadock files
	foreach my $mdock (@{$mdocks{$receptor}}){

		## Grabbing ligand from filename
		$ligand_count++;
		my ($ligand) = $mdock =~ /vs\.(\S+).out$/;

		if ($verbose){

			my $array_size = scalar @{$mdocks{$receptor}};
			if ($ligand_count == 1){
				print "  $array_size (total) receptor + ligand pairs for $receptor\n";
			}

			my $modulo = $ligand_count % 100;
			if ($modulo == 0){
				my $pad = length($array_size);
				$ligand_count = sprintf("%0${pad}d", $ligand_count);
				print "  Pair # $ligand_count of $array_size: $receptor (receptor) + $ligand (ligand)\n";
			}

		}

		## Running PPIscore + adding value to rankings
		unless (exists $rankings{$receptor}{$ligand}){
			my $ppi = `echo \$(ppiscore $mdock $npred)`;
			my ($ppiscore) = $ppi =~ /E = (\d+\.\d+)/;
			$rankings{$receptor}{$ligand} = $ppiscore;
		}
	}

	print TSV "\# $receptor\n";

	## Print PPIscores from best to worst for each receptor (easier to read)
	foreach my $ligand (reverse (sort { $rankings{$receptor}{$a} <=> $rankings{$receptor}{$b} } keys %{$rankings{$receptor}}) ){
		my $value = $rankings{$receptor}{$ligand};
		if ($value >= $ppicutoff){ 
			print TSV $receptor."\t".$ligand."\t".$value;
			if ($products){
				my ($locus) = $ligand =~ /(.*)-m\d$/;
				if (exists $products{$locus}){
					print TSV "\t".$products{$locus};
				}
				else {
					print TSV "\t".'N/A';
				}
			}
			print TSV "\n";
		}
	}

}
