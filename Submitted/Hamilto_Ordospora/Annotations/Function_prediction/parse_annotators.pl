#!/usr/bin/perl
## Pombert Lab, IIT, 2017
my $name = 'parse_annotators.pl';
my $version = '1.0';

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	This script parses the output of annotators to help assign putative functions to predicted proteins.
		Annotators are BLASTP searches against SwissProt/TREMBL and InterProScan 5 queries

USAGE		parse_annotators.pl -q BEOM2.proteins.queries \\
		-sl sprot.list -sb BEOM2.sprot.blastp.6 \\
		-tl trembl.list -tb BEOM2.trembl.blastp.6 \\
		-ip BEOM2.interpro.tsv

OPTIONS:
-q	List of proteins queried against annotators
-sl	List of proteins and their products in SwissProt
-sb	BLAST output file in outfmt 6 format of proteins queried against SwissProt
-tl	List of proteins and their products in TREMBL
-tb	BLAST output file in outfmt 6 format of proteins queried against TREMBL
-ip	TSV output from InterProScan

NOTE: The trembl.list file is large and will eat up at least 5 Gb of RAM
OPTIONS
die "$usage\n" unless @ARGV;

my $queries;
my $splist; my $spblast;
my $tblist; my $tbblast;
my $ipro;
GetOptions(
	'q=s' => \$queries,
	'sl=s' => \$splist,
	'sb=s' => \$spblast,
	'tl=s' => \$tblist,
	'tb=s' => \$tbblast,
	'ip=s' => \$ipro,
);

## Parsing SwissProt and TREMBL product lists
my $time = localtime(); my $tstart = time;
print "$time: Parsing SwissProt and TREMBL product lists $splist and $tblist, this might take a while...\n";
my %sprot; my %trembl;
open SP, "<$splist"; while (my $line = <SP>){chomp $line; my @cols = split("\t", $line); $sprot{$cols[0]}=$cols[1];}
open TB, "<$tblist"; while (my $line = <TB>){chomp $line; my @cols = split("\t", $line); $trembl{$cols[0]}=$cols[1];}
my $time_taken = time - $tstart; $time = localtime(); 
print "$time: Finished parsing SwissProt and TREMBL product lists $splist and $tblist in $time_taken seconds.\n";
close SP; close TB;

## Parsing SwissProt and TREMBL blast outfmt 6
$time = localtime(); 
print "$time: Parsing BLAST files $spblast and $tbblast...\n";
open SB, "<$spblast"; my %sphits;
while (my $line = <SB>){
	chomp $line;
	my @cols = split("\t", $line);	my $query =$cols[0]; my $hit = $cols[1]; my $evalue = $cols[10];
	if (exists $sphits{$query}){next;}
	elsif ($sprot{$hit} =~ /uncharacterized/i){next;} ## Discarding uninformative BLAST hists
	elsif ($sprot{$hit} =~ /hypothetical/i){next;} ## Discarding uninformative BLAST hists
	elsif ($sprot{$hit} =~ /predicted protein/i){next;} ## Discarding uninformative BLAST hists
	else{$sphits{$query}[0] = $sprot{$hit}; $sphits{$query}[1] = $evalue;}
}
close SB;
open TBB, "<$tbblast"; my %tbhits;
while (my $line = <TBB>){
	chomp $line;
	my @cols = split("\t", $line);	my $query =$cols[0]; my $hit = $cols[1]; my $evalue = $cols[10];
	if (exists $tbhits{$query}){next;}
	elsif ($trembl{$hit} =~ /uncharacterized/i){next;} ## Discarding uninformative BLAST hists
	elsif ($trembl{$hit} =~ /hypothetical/i){next;} ## Discarding uninformative BLAST hists
	elsif ($trembl{$hit} =~ /predicted protein/i){next;} ## Discarding uninformative BLAST hists
	else{$tbhits{$query}[0] = $trembl{$hit}; $tbhits{$query}[1] = $evalue;}
}
close TBB;
$time = localtime(); 
print "$time: Finished parsing BLAST files $spblast and $tbblast...\n";

## Parsing InterProScan 5 output
$time = localtime(); 
print "$time: Parsing InterProScan5 $ipro...\n";
open IP, "<$ipro";
my %pfam = (); my %hamap = ();
my %tigr = (); my %cdd = ();
while (my $line = <IP>){
	chomp $line;
	my @cols = split("\t", $line); ## Columns info from https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats
	## $cols[0] Protein Accession (e.g. P51587)
	## $cols[1] Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
	## $cols[2] Sequence Length (e.g. 3418)
	## $cols[3] Analysis (e.g. Pfam / PRINTS / Gene3D)
	## $cols[4] Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
	## $cols[5] Signature Description (e.g. BRCA2 repeat profile)
	## $cols[6] Start location
	## $cols[7] Stop location
	## $cols[8] Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
	## $cols[9] Status - is the status of the match (T: true)
	## $cols[10] Date - is the date of the run
	## $cols[11] (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprlookup option is switched on)
	## $cols[12] (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprlookup option is switched on)
	## $cols[13] (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
	## $cols[14] (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
	if ($cols[3] eq 'Pfam'){
		if (exists $pfam{$cols[0]}){if ($cols[8] < $pfam{$cols[0]}[1]){$pfam{$cols[0]}[0] = $cols[5]; $pfam{$cols[0]}[1] = $cols[8];}}
		else{$pfam{$cols[0]}[0] = $cols[5]; $pfam{$cols[0]}[1] = $cols[8];}
	}
	elsif ($cols[3] eq 'TIGRFAM'){
		if (exists $tigr{$cols[0]}){if ($cols[8] < $tigr{$cols[0]}[1]){$tigr{$cols[0]}[0] = $cols[5]; $tigr{$cols[0]}[1] = $cols[8];}}
		else{$tigr{$cols[0]}[0] = $cols[5]; $tigr{$cols[0]}[1] = $cols[8];}
	}
	elsif ($cols[3] eq 'Hamap'){
		if (exists $hamap{$cols[0]}){if ($cols[8] > $hamap{$cols[0]}[1]){$hamap{$cols[0]}[0] = $cols[5]; $hamap{$cols[0]}[1] = $cols[8];}}
		else{$hamap{$cols[0]}[0] = $cols[5]; $hamap{$cols[0]}[1] = $cols[8];}
	}
	elsif ($cols[3] eq 'CDD'){
		if (exists $cdd{$cols[0]}){if ($cols[8] > $cdd{$cols[0]}[1]){$cdd{$cols[0]}[0] = $cols[5]; $cdd{$cols[0]}[1] = $cols[8];}}
		else{$cdd{$cols[0]}[0] = $cols[5]; $cdd{$cols[0]}[1] = $cols[8];}
	}
}
close IP;
$time = localtime(); 
print "$time: Finished parsing InterProScan5 $ipro...\n";

## Creating the parsed annotations file
open QUE, "<$queries";
$queries =~ s/.proteins.queries$//;
$time = localtime(); 
print "$time: Writing annotations to $queries.annotations...\n"; $tstart = time;
open OUT, ">$queries.annotations";
while (my $line = <QUE>){
	chomp $line;
	print OUT "$line\t";
	if (exists $sphits{$line}){print OUT "$sphits{$line}[1]\t$sphits{$line}[0]\t";} else {print OUT "NA\thypothetical protein\t";}
	if (exists $tbhits{$line}){print OUT "$tbhits{$line}[1]\t$tbhits{$line}[0]\t";} else {print OUT "NA\thypothetical protein\t";}
	if (exists $pfam{$line}){print OUT "$pfam{$line}[1]\t$pfam{$line}[0]\t";} else {print OUT "NA\thypothetical protein\t";}
	if (exists $tigr{$line}){print OUT "$tigr{$line}[1]\t$tigr{$line}[0]\t";} else {print OUT "NA\thypothetical protein\t";}
	if (exists $hamap{$line}){print OUT "$hamap{$line}[1]\t$hamap{$line}[0]\t";} else {print OUT "NA\thypothetical protein\t";}
	if (exists $cdd{$line}){print OUT "$cdd{$line}[1]\t$cdd{$line}[0]\n";} else {print OUT "NA\tno motif found\n";}
}
close QUE, close OUT;
$time_taken = time - $tstart; $time = localtime(); 
print "$time: Done writing annotations in $time_taken seconds. Exiting...\n";
exit;
