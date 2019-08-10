#!/usr/bin/perl
## Pombert Lab, 2018
my $version = 0.2;
my $name = 'make_orthogroup_datasets.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Create Fasta datasets from Orthofinder Orthgroups output files.
USAGE		make_orthogroup_datasets.pl -f *.fasta -t *.tsv

-f (--fasta)	Multifasta files (proteins or nucleotides)
-t (--tsv)	Orthogroups.tsv file(s) from Orthofinder

NOTE: This script will split single- and multi-copy orthologs in distinct subfolders
OPTIONS
die "$usage\n" unless @ARGV;
my @fasta; my @csv;
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	't|tsv=s@{1,}' => \@csv
);

my %db; my $locus; ## Creating DB of sequences
while (my $fasta = shift@fasta){ 
	open FASTA, "<$fasta";
	$fasta =~ s/\.\w+$//; 
	while (my $line = <FASTA>){
		chomp $line;
		if ($line =~ /^>(\S+)/){
			$locus = $1; ## Keeping only the first part of the header
			$db{$fasta}{$locus}[0] = $fasta; ## Must keep track of file names in case fasta identifiers are identical between datatest and prevent overwrite of data...
			#print "$fasta\t$locus\n";
		} 
		else{$db{$fasta}{$locus}[1] .= $line;}
	}
}
my @OG; my @species;
while (my $csv = shift@csv){ ## Creating multifasta files for each ortthogroup
	print "\n"; system "dos2unix $csv"; print "\n"; ## Removing DOS format SNAFUs...
	open CSV, "<$csv"; open ORT, ">single_copy_$csv"; open PAR, ">multi_copy_$csv";
	$csv =~ s/.csv$//; $csv =~ s/.tsv$//;
	mkdir ("SINGLE_COPY_$csv",0755); mkdir ("MULTI_COPY_$csv",0755);
	my $mc = 0; my $sc = 0;
	while (my $line = <CSV>){
		#chomp $line; Chomping breaks the last tab with split(); -> chomping later instead
		if ($line =~ /^\t/){
			print ORT "Renamed_to\t$line";
			print PAR "Renamed_to\t$line";
			@species = split("\t", $line);
			for (1..$#species){print "Species # $_ = $species[$_]\n";}
		}
		else{
			@OG = split(/\t/, $line); ## @OG -> list of orthologs, single or multicopy
			my $og = $OG[0]; print "Working on $og...\n";
			if ($line =~ /,/){ ## Multicopy
				$mc++; $mc = sprintf("%05d", $mc);
				print PAR "MOG$mc\t$line\n";
				open OUT, ">./MULTI_COPY_$csv/MOG$mc.fasta"; ## SOG = Multi-copy orthogroups
				seq();
			}
			else { ## Single copy
				$sc++; $sc = sprintf("%05d", $sc);
				print ORT "SOG$sc\t$line\n";
				open OUT, ">./SINGLE_COPY_$csv/SOG$sc.fasta"; ## SOG = Single-copy orthogroups
				seq();
			}
		}
	}
}
sub seq{ ## Print sequence sub.
	for (1..$#species){
		my $organism = $species[$_]; chomp $organism;
		my $sp = $OG[$_];
		my @splits = split(",", $sp);
		while (my $para = shift@splits){
			chomp $para;
			if ((exists $db{$organism}{$para}[0])&&($db{$organism}{$para}[1] ne '')){
				print OUT ">$db{$organism}{$para}[0]\@$para\n";
				my @seq = unpack("(A60)*", $db{$organism}{$para}[1]);
				while (my $tmp = shift@seq){print OUT "$tmp\n";}
			}
		}
	}
	close OUT;
}
