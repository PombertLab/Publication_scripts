#!/usr/bin/perl
## Pombert Lab, IIT, 2017
my $name = 'curate_annotations.pl';
my $version = '1.0';

use strict; use warnings; use Getopt::Long qw(GetOptions); use FileHandle;

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Displays lists of functions predicted per proteins. User can select or enter desired annotation.
		Creates a tab-delimited .curated list of annotations.
USAGE		curate_annotations.pl -r -i BEOM2.annotations

OPTIONS:
-r	Resumes annotation from last curated locus_tag
-i	Input file (generated from parse_annotators.pl)
OPTIONS
die "$usage\n" unless @ARGV;

my $input;
my $resume;
GetOptions(
	'i=s' => \$input,
	'r' => \$resume
);

my %curated;
if ($resume){
	open RE, "<$input.curated";
	while (my $line = <RE>){
		chomp $line;
		if ($line =~ /^(\S+)\t(.*)$/){$curated{$1} = $2;}
	}
}
open IN, "<$input";
open OUT, ">>$input.curated"; OUT->autoflush(1); 

my $call; my $count = 0;
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^#/){next;}
	else{
		my @col = split("\t", $line); $count++; $count = sprintf("%04d", $count);
		## $col[0] => queries
		## $col[1] and $col[2] => SwissProt Evalues and Annotations
		## $col[3] and $col[4] => TREMBL Evalues and Annotations
		## $col[5] and $col[6] => Pfam Evalues and Annotations
		## $col[7] and $col[8] => TIGRFAM Evalues and Annotations
		## $col[9] and $col[10] => HAMAP Scores and Annotations
		## $col[11] and $col[12] => CDD Evalues and Motifs
		if (exists $curated{$col[0] }){next;}
		elsif (($col[2] eq 'hypothetical protein')&&($col[4] eq 'hypothetical protein')&&($col[6] eq 'hypothetical protein')&&($col[8] eq 'hypothetical protein')&&($col[10] eq 'hypothetical protein')&&($col[12] eq 'no motif found')){print OUT "$col[0]\thypothetical protein\n"; $curated{$col[0]} = 'in progress';}
		else{
			$curated{$col[0]} = 'in progress';
			print "\nPutative annotation(s) found for protein #$count: $col[0]:\n";
			print "1.\tSWISSPROT:\t$col[1]\t$col[2]\n";
			print "2.\tTREMBL:\t\t$col[3]\t$col[4]\n";
			print "3.\tPfam:\t\t$col[5]\t$col[6]\n";
			print "4.\tTIGRFAM:\t$col[7]\t$col[8]\n";
			print "5.\tHAMAP:\t\t$col[9]\t$col[10]\n";
			print "6.\tCDD:\t\t$col[11]\t$col[12]\n\n";
			print "Please enter selection [1-6] to assign annotation, [0] to annotate as 'hypothetical protein', [m] for manual annotation, or [x] to exit\n";
			chomp ($call = <STDIN>);
			if ($call eq 'x'){exit;}
			elsif ($call eq 'm'){
				print "Enter desired annotation: ";
				chomp (my $annot = <STDIN>);
				print OUT "$col[0]\t$annot\n";
			}
			elsif ($call =~ /^[0123456]$/){
				if ($call == 0){print OUT "$col[0]\thypothetical protein\n";}
				else{my $num = $call*2; print OUT "$col[0]\t$col[$num]\n";}
			}
			else {print "Wrong Input. Exiting to prevent SNAFUs...\n"; exit;}
		}
	}
}
