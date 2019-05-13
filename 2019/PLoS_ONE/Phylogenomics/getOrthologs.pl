#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2017

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = <<'OPTIONS';

Example: getOrthologs.pl -t Streptococcus_suis -p *.faa -o OrthologousGroups.csv -d output_dir

-t (--taxa)		## Desired name for the species/OTUs
-p (--proteins)		## Fasta file(s) of proteins 
-o (--ortho)		## Otholog CSV file created by OrthoFinder [Note: make sure that the file is in UNIX format]
-d (--dir)		## Output directory

OPTIONS
die "$usage" unless @ARGV;

## Defining variables
my $taxa;
my @fasta;
my $csv;
my $dir;
GetOptions(
	't|taxa=s' => \$taxa,
	'p|proteins=s@{1,}' => \@fasta,
	'o|ortho=s' => \$csv,
	'd|dir=s' => \$dir
);

## Creating hash of sequences; protein names must be unique
my %fasta = ();
my $key = undef;
while (my $fasta = shift@fasta){
	open FASTA, "<$fasta";
	while (my $line = <FASTA>){
		chomp $line;
		if ($line =~ /^>(\S+)/){$key = $1;}
		else {$fasta{$key} .= $line;}
	}
}

## Parsing Orthofinder output/creating fasta files
open CSV, "<$csv";
my @species = ();
while (my $line = <CSV>){
	if ($line =~ "\t\n"){next;}		## Looking for missing protein at end of line; if missing skip to next line
	elsif ($line =~/,/){next;} 		## Skipping orthologous group when paralogs are found
	elsif ($line =~/\t\t/){next;}	## Skipping blank orthologs
	else {
		chomp $line;
		if ($line =~ /^\t/){$line =~ s/.faa//g; $line =~ s/.fasta//g; @species = split("\t", $line);}	## Grabbing species names
		elsif ( $line =~ /^OG/){
			my @OG = split("\t", $line); ## OG = $OG[0];
			open OUT, ">$OG[0].orth.fasta";
			my @sequence = ();
			for (1..$#OG){
				print OUT ">${taxa}_$species[$_]\@$OG[$_]\n";
				@sequence = unpack("(A80)*", $fasta{$OG[$_]});
				while (my $seq = shift@sequence){print OUT "$seq\n";}
			}
		}
	}
}
close CSV; close OUT;

## Cleaning up
if ($dir){system "mkdir $dir; mv *.orth.fasta $dir/";}