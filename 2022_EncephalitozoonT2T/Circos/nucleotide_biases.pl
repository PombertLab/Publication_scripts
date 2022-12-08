#!/usr/bin/perl
## Pombert Lab, 2022
my $name = 'nucleotide_biases.pl';
my $version = '0.3b';
my $updated = '2022-06-26';

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
		 distributions for easy plotting with MS Excel or other tool. Can also generate
		 coordinate files for Circos.

COMMAND		$name \\
		  -f *.fasta \\
		  -o output_directory \\
		  -w 1000 \\
		  -s 500 \\
		  -c

-f (--fasta)	Fasta file(s) to process
-o (--outdir)	Output directory [Default: ntBiases]
-w (--winsize)	Sliding window size [Default: 1000]
-s (--step)		Sliding window step [Default: 500]
-c (--circos)	Output files for Circos plotting
OPTIONS
die "\n$usage\n" unless @ARGV;

my @fasta;
my $outdir = 'ntBiases';
my $winsize = 1000;
my $step = 500;
my $circos;
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
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

### Hash to store percentage values
my %percent;

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

	my $fasta_dir = "$outdir/$fileprefix";
	unless (-d $fasta_dir){
		mkdir ($fasta_dir,0755) or die "Can't create $fasta_dir: $!\n";
	}

	### Circos
	no warnings 'once'; ## Filehandles are used in a loop...
	my @filehandles = (*GC, *AT, *GT, *AC, *GA, *CT);

	my $circos_GC;
	my $circos_AT;
	my $circos_GT;
	my $circos_AC;
	my $circos_GA;
	my $circos_CT;
	
	my $circos_kar;
	my $id = 0;

	if ($circos){

		my $circos_dir = $outdir.'/CIRCOS/'.$fileprefix;
		unless (-d $circos_dir) {
				make_path($circos_dir,{mode => 0755}) or die "Can't create $circos_dir: $!\n";
		}

		## Creating a "karyotype" file for Circos
		$circos_kar =   $circos_dir.'/'.$fileprefix.'.genotype';
		open KAR, ">", $circos_kar or die "Can't create $circos_kar: $!\n";
		print KAR '#chr - ID LABEL START END COLOR'."\n";
		
		## Nucleotide biases
		for my $fh (@filehandles){
			my ($lfh) = $fh =~ /(\w+)$/; ## grabbing test from filehandle
			my $filename = $circos_dir.'/'.$fileprefix.'.'.$lfh;
			open $fh, ">", $filename or die "Can't create $filename: $!\n";
			print $fh '#chr START END GC_PERCENTAGE'."\n";
		}

	}

	### Iterating through each sequence in the FASTA file
	foreach my $sequence (sort (keys %sequences)){

		my $outfile = $fasta_dir.'/'.$fileprefix.'.'.$sequence.'.tsv';
		open BIAS, ">", $outfile or die "Can't create $outfile: $!\n";

		print BIAS "# Location\t% GC\t% AT\t% AG\t% CT\t% GT\t% AC\n";

		### Sliding windows
		my $seq = $sequences{$sequence};
		my $csize = length $seq;

		if ($circos){
			my $terminus = $csize - 1;
			$id++;
			print KAR "chr - $sequence $id 0 $terminus black\n";
		}

		my $x;
		for ($x = 0; $x <= ($csize - $winsize); $x += $step){

			my $subseq = substr($seq, $x, $winsize);
			my $end = $x + $winsize - 1;
			
			%percent = ();
			biases($subseq,$x);

			if ($circos){
				foreach my $fh (@filehandles){
					my ($lfh) = $fh =~ /(\w+)$/;
					print $fh "$sequence $x $end $percent{$lfh}\n";
				}
			}

		}

		### Working on leftover string < $winsize
		my $modulo = $csize % $winsize;
		my $subseqleft = substr($seq, -$modulo, $modulo);
		my $leftover_size = length $subseqleft;

		%percent = ();
		biases($subseqleft,$x);

		if ($circos){
			my $end = $csize - 1;
			foreach my $fh (@filehandles){
				my ($lfh) = $fh =~ /(\w+)$/;
				print $fh "$sequence $x $end $percent{$lfh}\n";
			}
		}

		close BIAS;

	}

	close FASTA;

}

### Subroutine(s)

sub biases {

		my $curseq = $_[0];
		my $pos = $_[1];

		my $gc = $curseq =~ tr/GgCc//;
		my $at = $curseq =~ tr/AaTt//;
		my $ga = $curseq =~ tr/GgAa//;
		my $ct = $curseq =~ tr/CcTt//;
		my $gt = $curseq =~ tr/GgTt//;
		my $ac = $curseq =~ tr/AaCc//;
		

		$percent{'GC'} = $gc = ($gc/$winsize) * 100;
		$percent{'AT'} = $at = ($at/$winsize) * 100;
		$percent{'GA'} = $ga = ($ga/$winsize) * 100;
		$percent{'CT'} = $ct = ($ct/$winsize) * 100;
		$percent{'GT'} = $gt = ($gt/$winsize) * 100;
		$percent{'AC'} = $ac = ($ac/$winsize) * 100;

		print BIAS "$pos\t$gc\t$at\t$ga\t$ct\t$gt\t$ac\n";

}