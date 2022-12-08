#!/usr/bin/perl
## Pombert Lab, 2022
my $name = 'nucleotide_biases.pl';
my $version = '0.2';
my $updated = '2022-06-22';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $usage = <<"OPTIONS";
NAME		$name
VERSION		$version
UPDATED		$updated
SYNOPSIS	Generates tab-delimited sliding windows of GC, AT, purine and pyrimidine
		 distributions for easy plotting with MS Excel or other tool.

COMMAND		$name \\
		  -f *.fasta \\
		  -o output_directory \\
		  -w 1000 \\
		  -s 500

-f (--fasta)	Fasta file(s) to process
-o (--outdir)	Output directory [Default: ntBiases]
-w (--winsize)	Sliding window size [Default: 1000]
-s (--step)		Sliding window step [Default: 500]
OPTIONS
die "\n$usage\n" unless @ARGV;

my @fasta;
my $outdir = 'ntBiases';
my $winsize = 1000;
my $step = 500;
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'o|outdir=s' => \$outdir,
	'w|winsize=i' => \$winsize,
	's|step=i' => \$step
);

### Check if output directory / subdirs can be created
$outdir =~ s/\/$//;
unless (-d $outdir) {
	make_path( $outdir, { mode => 0755 } )  or die "Can't create $outdir: $!\n";
}

### Iterating through FASTA file(s)
while (my $fasta = shift@fasta){

	my ($basename) = fileparse($fasta);
	my ($fileprefix) = $basename =~ /(\S+)\.\w+$/;

	open FASTA, "<", $fasta or die "Can't open $fasta: $!\n";

	### Creating database of sequences (could be multifasta)
	my %sequences;
	my $seqname;

	while (my $line = <FASTA>){
		chomp $line;
		if ($line =~ /^>(\S+)/){
			$seqname = $1;
		}
		else {
			$sequences{$seqname} .= $line;
		}
	}

	### Iterating through each sequence in the FASTA file
	foreach my $sequence (sort (keys %sequences)){

		my $outfile = $outdir.'/'.$fileprefix.'.'.$sequence.'.tsv';
		open GC, ">", $outfile or die "Can't create $outfile: $!\n";

		print GC "# Location\t% GC\t% AT\t% Purines\t% Pyrimidines\t% GT\t% AC\n";

		### Sliding windows
		my $seq = $sequences{$sequence};
		my $csize = length $seq;
		my $x;
		for ($x = 0; $x <= ($csize - $step); $x += $step){
			my $subseq = substr($seq, $x, $winsize);
			my $gc = $subseq =~ tr/GgCc//;
			my $at = $subseq =~ tr/AaTt//;
			my $pur = $subseq =~ tr/GgAa//;
			my $pyr = $subseq =~ tr/CcTt//;
			my $gt = $subseq =~ tr/GgTt//;
			my $ac = $subseq =~ tr/AaCc//;
			$gc = ($gc)/$winsize * 100;
			$at = ($at)/$winsize * 100;
			$pur = ($pur)/$winsize * 100;
			$pyr = ($pyr)/$winsize * 100;
			$gt = ($gt)/$winsize * 100;
			$ac = ($ac)/$winsize * 100;
			print GC "$x\t$gc\t$at\t$pur\t$pyr\t$gt\t$ac\n";
		}

		### Working on leftover string < $winsize
		my $modulo = $csize % $winsize;
		my $subseqleft = substr($seq, -$modulo, $modulo);
		my $leftover_size = length $subseqleft;
		my $gc = $subseqleft =~ tr/GgCc//;
		my $at = $subseqleft =~ tr/AaTt//;
		my $pur = $subseqleft =~ tr/GgAa//;
		my $pyr = $subseqleft =~ tr/CcTt//;
		my $gt = $subseqleft =~ tr/GgTt//;
		my $ac = $subseqleft =~ tr/AaCc//;
		$gc = ($gc)/$leftover_size * 100;
		$at = ($at)/$winsize * 100;
		$pur = ($pur)/$winsize * 100;
		$pyr = ($pyr)/$winsize * 100;
		$gt = ($gt)/$winsize * 100;
		$ac = ($ac)/$winsize * 100;
		print GC "$x\t$gc\t$at\t$pur\t$pyr\t$gt\t$ac\n";

		close GC;

	}

	close FASTA;

}