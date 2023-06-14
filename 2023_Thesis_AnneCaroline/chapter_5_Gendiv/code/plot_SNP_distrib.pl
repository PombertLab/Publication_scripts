#!/usr/bin/perl
## Pombert Lab, 2023
my $name = 'plot_SNP_distrib.pl';
my $version = '0.1a';
my $updated = '2023-04-26';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Cwd qw(abs_path);

my $options = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Plots variants distribution per contig from VCF files with Circos (using sliding windows)

COMMAND		${name} \\
		  -v *.vcf \\
		  -f *.fasta \\
		  -o outdir/

OPTIONS:
-f (--fasta)	Input files in fasta format
-o (--outdir)	Output directory [Default: ./]
-v (--vcf)	Input files in vcf format
-s (--step)	Size of the steps between windows [Default: 500]
-w (--window)	Width of the sliding windows [Default: 1000]
OPTIONS
die "\n$options\n" unless @ARGV;

my @vcf;
my @fasta;
my $outdir = './';
my $slide = 500;
my $window = 1000;

GetOptions(
	'v|vcf=s@{1,}' => \@vcf,
	'f|fasta=s@{1,}' => \@fasta,
	'o|outdir=s' => \$outdir,
	's|slide=i' => \$slide,
	'w|window=i' => \$window
);

die "\nPlease enter at least one fasta file...\n\n" unless @fasta;

## Checking outdir
$outdir = abs_path($outdir);
unless (-d $outdir){
	mkdir ($outdir, 0755) or die $!;
	mkdir ("$outdir/data", 0755) or die $!;
	mkdir ("$outdir/images", 0755) or die $!;
}

## Populating list of SNPs
my %db;
my $tot_snps = 0;

while (my $vcf = shift @vcf){

	open IN, '<', $vcf or die "Can\'t open $vcf: $!\n";

	while (my $line = <IN>){
		
		chomp $line;

		if ($line =~ /^#/){
			next;
		}

		elsif ($line =~ /^(\S+)\s+(\d+).*HET=(\d);HOM=(\d)/){

			my $locus = $1;
			my $position = $2;
			my $het = $3;
			my $hom = $4;

			if ($het != 0){
				$db{$locus}{'snp'}{$position} = 'het';
				$tot_snps ++;
			}

			elsif ($hom != 0){
				$db{$locus}{'snp'}{$position} = 'hom';
				$tot_snps ++;
			}
		}
	}
	close IN;
}

## Populating list of sequences
my @seq_list;
my $seq;

while (my $fasta = shift @fasta) {

	open IN, '<', $fasta or die "Can\'t open $fasta: $!\n";

	while (my $line = <IN>){

		chomp $line;

		if ($line =~ /^>(\S+).*$/){
			$seq = $1;
			push (@seq_list, $seq);
		}

		else {
			$db{$seq}{'seq'} .= $line;
		}

	}
	close IN;
}


my $mmax;
my $mmin;

## Ideogram for Circos
my $ideogram = $outdir.'/data/ideogram.conf';
open my $id, '>', $ideogram or die $!;

my $ideogram_data = <<'IDEO';
<ideogram>

<spacing>
default = 0.01r
</spacing>

radius    = 0.9r
thickness = 5p
fill      = yes
show_label       = yes
label_font       = bold 
label_radius     = dims(ideogram,radius) + 0.07r
label_size       = 36
label_parallel   = yes

</ideogram>
IDEO
print $id $ideogram_data."\n";
close $id;

############## Ticks
my $ticks = $outdir.'/data/ticks.conf';
open my $tk, '>', $ticks or die $!;

my $ticks_data = <<'TICKS';
show_ticks          = yes
show_tick_labels    = yes

<ticks>
radius           = 1r
orientation	= out
tick_separation  = 3p
label_separation = 1p
multiplier       = 1e-6
color            = black
thickness        = 3p

<tick>
size             = 20p
spacing        = 10u
show_label     = yes
label_size     = 15p
label_offset   = 5p
format         = %d
suffix         = kb
grid           = yes
grid_color     = dgrey
grid_thickness = 1p
grid_start     = 1r
grid_end       = 1.1r
</tick>

<tick>
size             = 10p
spacing        = 5u
show_label     = yes
label_size     = 10p
label_offset   = 5p
format         = %d
thickness      = 1.5p
color          = vdgrey
</tick>

<tick>
size             = 5p
spacing        = 1u
show_label     = no
thickness      = 1p
color          = dgrey
</tick>

</ticks>
TICKS
print $tk $ticks_data."\n";
close $tk;

for my $contig (keys %db){

	my $karyotype = $outdir.'/data/'.$contig.'.karyotype';
	my $config =  $outdir.'/data/'.$contig.'.conf';
	my $snpfile = $outdir.'/data/'.$contig.'.snp';
	
	my $pngdir = $outdir.'/images/';
	my $image = $contig.'.png';

	## Creating a karyotype for Circos
	open my $ka, '>', $karyotype or die $!;
	print $ka '#chr - ID LABEL START END COLOR'."\n";
	my $seqlen = length($db{$contig}{'seq'}) - 1;
	print $ka "chr - $contig $contig 0 $seqlen grey"."\n";
	print $ka "chr - spacer sp 0 100 grey"."\n";

	## Creating a conf file for Circos
	open my $cg, '>', $config  or die $!;

	print $cg "<<include $ideogram>>"."\n";
	print $cg "<<include $ticks>>"."\n\n";
	print $cg "karyotype = $karyotype"."\n";
	print $cg 'chromosomes_units=10000'."\n\n";

	print $cg '<plots>'."\n";
	print $cg '<plot>'."\n";
	print $cg 'type      = heatmap'."\n";
	print $cg "file    = $snpfile"."\n";
	print $cg 'color   = purples-13'."\n";
	print $cg 'min     = 0'."\n";
	print $cg 'max     = 30'."\n";
	print $cg 'r1      = 0.99r'."\n";
	print $cg 'r0      = 0.29r'."\n";
	print $cg 'log_scale_base	= 2'."\n";
	print $cg '</plot>'."\n";
	print $cg '</plots>'."\n\n";

	print $cg '<image>'."\n";
	print $cg '<<include etc/image.conf>>'."\n";
	print $cg '</image>'."\n\n";
	print $cg '<<include etc/colors_fonts_patterns.conf>>'."\n";
	print $cg '<<include etc/housekeeping.conf>>'."\n";

	## Creating a "SNP" file for Circos
	open my $fh, '>', $snpfile or die $!;
	print $fh '#chr START END SNP_count'."\n";

	## Working on sequences
	my $len = length($db{$contig}{'seq'});
	my $terminus = $len - 1;
	my $max;
	my $min;
	
	print "Working on $contig ... ";
	
	for (my $i = 0; $i <= $len - ($window-1); $i += $slide) {

		my $snp = 0;
		
		for ($i..($i+$window-1)){
			if (exists $db{$contig}{'snp'}{$_}){
				$snp++;
			}
		}
		
		my $end = $i + $window - 1;
		print $fh "$contig $i $end $snp\n";
		
		if (!defined $max){$max = $snp;}
		if (!defined $mmax){$mmax = $snp;}
		if (!defined $min){$min = $snp;}
		if (!defined $mmin){$mmin = $snp;}
		if ($snp > $max){$max = $snp;}
		if ($snp > $mmax){$mmax = $snp;}
		if ($snp < $min){$min = $snp;}
		if ($snp < $mmin){$mmin = $snp;}
	
	}

	print "Min/Max SNP per window = $min / $max\n";
	
	close $ka;
	close $cg;
	close $fh;

	#### Running circos

	system "circos \\
	  -conf $config \\
	  -outputfile $image \\
	  -outputdir $pngdir";

}

print "\nMin SNP per window = $mmin\n";
print "Max SNP per window = $mmax\n";
print "Total number of SNPs = $tot_snps\n";

exit;
