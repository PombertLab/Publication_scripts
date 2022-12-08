#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2022
my $name = 'pombase_matches.pl';
my $version = '0.2';
my $updated = '2022-07-30';

use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(GetOptions);

my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Parses the output from pombase BLAST/FOLDSEEK searches against 
		target organism

OPTIONS:
-f (--fseek)	Foldseek output file(s)
-b (--blast)	Blast output file (in oufmt 6 format)
-e (--evalue)	E-value cutoff [Default: 1e-05]
-l (--list)	Tab-delimited list of GOs and related genes to use as scaffold for output table
-p (--prod)	Ta-delimited list of locus tags and their products for the target organism
-o (--out)	Output table in TSV format
-q (--qdir)		Folder containing PDB files used as query
-s (--sdir)		Folder containing PDB files used as subject
-m (--mican)	Folder to store Mican results

USAGE
die "\n$usage\n" unless @ARGV;

my @fseek;
my $blast;
my $eval_cutoff = '1e-05';
my $list;
my $prod;
my $outfile;
my $qdir = './Data/PDB';
my $sdir = './Data/ALPHAFOLD_3D_Parsed';
my $mican_dir = './Mican';
GetOptions(
	'f|fseek=s@{1,}' => \@fseek,
	'b|blast=s' => \$blast,
	'e|evalue=s' => \$eval_cutoff,
	'l|list=s' => \$list,
	'p|prod=s' => \$prod,
	'o|out=s' => \$outfile,
	'q|qdir=s' => \$qdir,
	's|sdir=s' => \$sdir,
	'm|mican=s' => \$mican_dir
);

my %hits;
my %products;

## Mican dir
unless (-d $mican_dir){
	mkdir ($mican_dir,0755) or die "Can't create $mican_dir: $!\n";
}

## Product list
if ($prod){
	open PROD, "<", $prod or die "Can't open $prod: $!\n";
	while (my $line = <PROD>){
		chomp $line;
		my ($locus, $product) = $line =~ /^(\w+)\t(.*)$/;
		$products{$locus} = $product;
	}
}

## Foldseek results
if (@fseek){

	while (my $fseek = shift@fseek){

		my $basename = fileparse($fseek);
		my ($gene) = $basename =~ /^([^_]+)/;
		$gene = lc($gene);

		open FSEEK, "<", $fseek or die "Can't open $fseek: $!\n";
		while (my $line = <FSEEK>){

			chomp $line;

			unless ($line =~ /^$/){

				my @data = split ("\t", $line);
				my $query = $data[0];
				my $subject = $data[1];
				my ($match) = $data[1] =~ /^(\w+)/;
				my $evalue = $data[-2];

				if ($evalue <= $eval_cutoff){
					if (exists $hits{$gene}{$match}{'fseek'}){
						next;
					}
					else {
						$hits{$gene}{$match}{'fseek'} = $evalue;

						## Mican
						my $mican_outfile = $mican_dir.'/'.$query.'.vs.'.$subject.'.mican';
						unless (-e $mican_outfile){
							my $mican = `mican -x $qdir/$query $sdir/$subject -a $mican_outfile`;
						}

						open MICAN, "<", $mican_outfile;
						while (my $line = <MICAN>){
							chomp $line;
							if ($line =~ /\# sTM-score\s+(\S+)/){
								my $sTMscore = $1;
								print "Mican sTM-score $query vs $subject = $sTMscore\n";
								$hits{$gene}{$match}{'tm'} = $sTMscore;
							}
						}
						close MICAN;
					}
				}

			}

		}
		close FSEEK;

	}

}

## Blastp results
if ($blast){

	open BLAST, "<", $blast or die "Can't open $blast: $!\n";
	while (my $line = <BLAST>){

		chomp $line;
		unless ($line =~ /^$/){

			my @data = split ("\t", $line);
			my ($gene) = $data[0] =~ /^([^_]+)/;
			$gene = lc($gene);
			my $match = $data[1];
			my $evalue = $data[-2];

			if ($evalue <= $eval_cutoff){
				if (exists $hits{$gene}{$match}{'blast'}){
					next;
				}
				else {
					$hits{$gene}{$match}{'blast'} = $evalue;
				}
			}

		}

	}
	close BLAST;

}

## Creating output table (tab-delimited)
open LIST, "<", $list or die "Can't open $list: $!\n";
open OUT, ">", $outfile or die "Can't create $outfile: $!\n";
print OUT "### Gene\tDescription\tMatch\tGenbank annotation\tBLAST evalue\tFoldseek evalue\tsTM-score\n\n";

while (my $line = <LIST>){

	chomp $line;
	my $gene;

	if ($line =~ /^#/){
		print OUT "$line\n";
	}
	elsif ($line =~ /^$/){
		print OUT "$line\n";
	}
	else {

		my ($gene, $desc) = split ("\t", $line);

		if (exists $hits{$gene}){

			foreach my $match (sort (keys %{$hits{$gene}})){

				print OUT "$gene\t$desc\t$match\t$products{$match}\t";

				if (exists $hits{$gene}{$match}{'blast'}){
					print OUT "$hits{$gene}{$match}{'blast'}\t";
				}
				else {
					print OUT "---\t";
				}
				if (exists $hits{$gene}{$match}{'fseek'}){
					my $evalue = $hits{$gene}{$match}{'fseek'};
					my $tmscore = $hits{$gene}{$match}{'tm'};
					print OUT "$evalue\t";
					print OUT "$tmscore\n";
				}
				else {
					print OUT "---\t---\n";
				}
			}

		}
		else {
			print OUT "$gene\t$desc\tno match\t---\t---\t---\t---\n";
		}
		print OUT "\n";
	}

}