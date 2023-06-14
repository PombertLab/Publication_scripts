#!/usr/bin/env perl
## Pombert Lab, Illinois Tech, 2013
my $name = 'create_CC_fasta.pl';
my $version = '0.1';
my $updated = '2023-05-05';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;

my $usage = <<"USAGE";

COMMAND	${name} -l cc_list.txt -f *.fasta -o outdir

USAGE
die "\n$usage\n" unless @ARGV;

my $cc_list = 'cc_list.txt';
my @fasta;
my $outdir = './';
GetOptions(
	'l|list=s' => \$cc_list,
	'f|fasta=s@{1,}' => \@fasta,
	'o|out=s' => \$outdir
);

## Output directory
unless (-d $outdir){
	mkdir ($outdir, 0755) or die $!;
}

### Creating a database of sequences
my %sequences;

while (my $fasta = shift @fasta){

	open FASTA, '<', $fasta or die $!;
	my ($basename,$path) = fileparse($fasta);
	$basename =~ s/\.\w+$//; ## Removing file extension

	my $seqname;

	while (my $line = <FASTA>){

		chomp $line;

		if ($line =~ /^>(\S+)/){
			$seqname = $1;
		}
		else {
			$sequences{$basename}{$seqname} .= $line;
		}

	}

}

### Creating a database of chromosome cores
my %cc;

open CC, '<', $cc_list or die $!;

while (my $line = <CC>){

	chomp $line;
	
	my @data = split("\t", $line);
	my $locus = $data[0];
	my $start = $data[1];
	my $end = $data[2];

	$cc{$locus}{'start'} = $start;
	$cc{$locus}{'end'} = $end;

}

### Creating CC fasta files

foreach my $fasta (keys %sequences){

	my $outfile = $outdir.'/'.$fasta.'.cc.fasta';
	open my $fh, '>', $outfile or die $!;

	foreach my $chromosome ( keys %{$sequences{$fasta}} ){

		my $start = $cc{$chromosome}{'start'} - 1;
		my $end = $cc{$chromosome}{'end'} - 1;
		my $cc_substring = substr($sequences{$fasta}{$chromosome}, $start, $end);

		print $fh '>'.$chromosome." \[CC subset: start: $start; end $end\]\n";
		my @seq = unpack ("(A60)*", $cc_substring);
		while (my $tmp = shift@seq){
			print $fh "$tmp\n";
		}

	}

}