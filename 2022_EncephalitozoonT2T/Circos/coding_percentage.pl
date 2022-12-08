#!/usr/bin/perl
## Pombert Lab, 2022
my $name = 'coding_percentage.pl';
my $version = '0.1';
my $updated = '2022-07-02';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $usage = <<"OPTIONS";
NAME		$name
VERSION		$version
UPDATED		$updated
SYNOPSIS	Generates tab-delimited sliding windows of coding percentages.
		 Useful to find loci of lower coding density. Can also generate
		 coordinate files for Circos.

COMMAND		$name \\
		  -f *.gb \\
		  -o output_directory \\
		  -w 1000 \\
		  -s 500 \\
		  -c

-g (--gbk)		GenBank file(s) to process
-o (--outdir)	Output directory [Default: cdPercent]
-w (--winsize)	Sliding window size [Default: 1000]
-s (--step)		Sliding window step [Default: 500]
-c (--circos)	Output files for Circos plotting
OPTIONS
die "\n$usage\n" unless @ARGV;

my @gbk;
my $outdir = 'cdPercent';
my $winsize = 1000;
my $step = 500;
my $circos;
GetOptions(
	'f|fasta=s@{1,}' => \@gbk,
	'o|outdir=s' => \$outdir,
	'w|winsize=i' => \$winsize,
	's|step=i' => \$step,
	'c|circos' => \$circos
);

### Check if output directory / subdirs can be created
$outdir =~ s/\/$//;
unless (-d $outdir) {
	make_path($outdir,{mode => 0755})  or die "Can't create $outdir: $!\n";
}

### Iterating through GBK file(s)
while (my $gbk = shift@gbk){

	my ($basename) = fileparse($gbk);
	my ($fileprefix) = $basename =~ /(\S+)\.\w+$/;

	open GBK, "<", $gbk or die "Can't open $gbk: $!\n";

	### Creating database of sequences (could be a multigbk)
	my %loci;
	my $locus;
	my $locus_size;

	while (my $line = <GBK>){

		chomp $line;
		
		if ($line =~ /^LOCUS\s+(\S+)\s+(\d+)\s+bp/){
		
			$locus = $1;
			$locus_size = $2;

			for (1..$locus_size){
				$loci{$locus}{$_} = 0;
			}

		}
		elsif ($line =~ /^^\s+gene\s+<?(\d+)\.\.>?(\d+)/){
			
			my $start = $1;
			my $end = $2;

			for ($start..$end){
				$loci{$locus}{$_} = 1;
			}

		}
		elsif ($line =~ /^\s+gene\s+complement\(<?(\d+)\.\.>?(\d+)/) {

			my $start = $1;
			my $end = $2;
	
			for ($start..$end){
				$loci{$locus}{$_} = 1;
			}
		}
	}

	### Circos
	if ($circos){
		my $circos_file = $outdir.'/'.$fileprefix.'.cdpercent.circos';
		open CIRCOS, ">", $circos_file or die "Can't create $circos: $!\n";
		print CIRCOS '#chr START END CODING_PERCENTAGE'."\n";
	}

	### Iterating through each sequence in the GBK file
	foreach my $sequence (keys %loci){

		my $len = scalar (keys %{$loci{$sequence}});
		my $outfile = $outdir.'/'.$fileprefix.'.'.$sequence.'.tsv';
		open TSV, ">", $outfile or die "Can't create $outfile: $!\n";
		print TSV "# Location\t% coding\n";

		### Sliding windows
		my $x;
		for ($x = 1; $x <= ($len - $winsize - 1); $x += $step){

			my $sum = 0;
			my $end = $x + $winsize - 1;
			for my $pos ($x..$end){
				$sum += $loci{$sequence}{$pos};
			}

			my $percent = ($sum/$winsize)*100;
			$percent = sprintf("%.1f", $percent);
			print TSV "$x\t$percent\n";

			if ($circos){
				my $cstart = $x - 1;
				my $cend = $end - 1;
				print CIRCOS "$sequence $cstart $cend $percent\n";
			}

		}

		### Working on leftover chunk
		my $leftover_size = $len - $x - 1;

		my $sum = 0;
		for my $pos ($x..$len){
			$sum += $loci{$sequence}{$pos};
		}

		my $percent = ($sum/$winsize)*100;
		$percent = sprintf("%.1f", $percent);
		print TSV "$x\t$percent\n";

		if ($circos){
			my $cstart = $x - 1;
			my $cend = $len - 1;
			print CIRCOS "$sequence $cstart $cend $percent\n";
		}

		close TSV;

	}

	close CIRCOS;

}
