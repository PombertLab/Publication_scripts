#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2022

my $name = 'check_for_telomeres.pl';
my $version = '0.4';
my $updated = '2022-10-14';

##### Pragmas + modules ########################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);

###### Command line options ####################################################

my $usage =<<"USAGE";
NAME		${name}
VERSION		${version}
UPPDATED	${updated}
SYNOPSIS	Checks the first and last X nucleotides of contigs for the presence 
		of telomere repeat units and calculates the top 5-mers and 6-mers

COMMAND		${name} \\
		  -f file.fasta \\
		  -o outdir \\
		  -w 1000 \\
		  -p TTAGG

OPTIONS:
-f (--fasta)	FASTA file to check
-o (--outdir)	Output dir location [Default: ./]
-w (--window)	Window size (first/last x nucleotides) to check [Default: 1000]
-p (--pattern)	Telomeric pattern to search for [Default: TTAGG]
-m (--min)	Minimum number of repeats [Default: 10]
USAGE
die "\n$usage\n" unless @ARGV;

my $fasta;
my $outdir = './';
my $window = 1000;
my $pattern = 'TTAGG';
my $min = 10;
GetOptions(
	'f|fasta=s' => \$fasta,
	'o|outdir=s' => \$outdir,
	'w|window=i' => \$window,
	'p|pattern=s' => \$pattern,
	'm|min=i' => \$min
);

###### I/O #####################################################################

open FASTA, "<", $fasta or die "Can't open $fasta: $!\n";

unless (-d $outdir){
	make_path($outdir, { verbose => 0, mode => 0755});
}
my $out_check = $outdir.'/'.'telomeres_check.txt';
my $out_summary = $outdir.'/'.'telomeres_summary.txt';
my $out_kmers = $outdir.'/'.'telomeres_kmers.txt';
open CHK, ">", $out_check or die "Can't create $out_check: $!\n";
open SUM, ">", $out_summary or die "Can't create $out_summary: $!\n";
open KMER, ">", $out_kmers or die "Can't create $out_kmers: $!\n";

##### Creating database of sequences ###########################################

my %sequences;
my $seqname;

while (my $line = <FASTA>){
	chomp $line;
	if ($line =~ /^>(.*)$/){
		$seqname = $1;
	}
	else {
		$sequences{$seqname} .= uc($line);
	}
}

##### Calculating kmers ########################################################

my %kmers;

foreach my $contig (keys %sequences){

	my $sequence = $sequences{$contig};
	my $length = length($sequence);

	my $boundary = $length - $window;

	### Kmers ranging from 5 to 8 should cover most possible telomere repeat units
	for my $klen (5..8){
		for (my $i = 0; $i <= ($length - $boundary - $klen); $i++){
			my $kmer = substr($sequence, $i, $klen);
			$kmers{$klen}{$kmer}++;
		}
		for (my $i = $boundary; $i <= ($length - $klen); $i++){
			my $kmer = substr($sequence, $i, $klen);
			$kmers{$klen}{$kmer}++;
		}
	}

}

for my $klen (5..8){

	my $counter = 0;
	print KMER "### Top $klen-mers\n";
	foreach my $key (reverse (sort { $kmers{$klen}{$a} <=> $kmers{$klen}{$b} } keys %{$kmers{$klen}})) {
		$counter++;
		print KMER $counter."\t".$key."\t".$kmers{$klen}{$key}."\n";
		last if ($counter == 10);
	}
	print KMER "\n";

}

##### Working on summary #######################################################

my $number_of_seq = scalar(keys %sequences);
my @lengths;
my $sum_len = 0;

foreach my $contig (keys %sequences){
	my $sequence = $sequences{$contig};
	my $length = length($sequence);
	$sum_len += $length;
	push(@lengths, $length);
}

@lengths = sort {$b <=> $a} @lengths;

# Average
my $longest = commify($lengths[0]);
my $shortest = commify($lengths[-1]);
my $average_len = commify(int($sum_len/$number_of_seq));

# median
my $median;
my $median_pos = $number_of_seq/2;
if ($median_pos =~ /^\d+$/){ $median = $lengths[$median_pos]; }
else {
	my $med1 = int($median_pos);
	my $med2 = $med1 + 1;
	$median = (($lengths[$med1] + $lengths[$med2])/2);
}
$median = sprintf("%.0f", $median);
$median = commify($median);

### N50, N75, N90
# thresholds to reach for N metrics
my $n50_td = $sum_len*0.5;
my $n75_td = $sum_len*0.75;
my $n90_td = $sum_len*0.9;
# sums to calculate
my $nsum50 = 0;
my $nsum75 = 0;
my $nsum90 = 0;
# metrics to capture
my $n50;
my $n75,
my $n90;

foreach (@lengths){
	$nsum50 += $_;
	if ($nsum50 >= $n50_td){
		$n50 = $_;
		last;
	}
}
foreach (@lengths){
	$nsum75 += $_;
	if ($nsum75 >= $n75_td){
		$n75 = $_;
		last
	}
}
foreach (@lengths){
	$nsum90 += $_;
	if ($nsum90 >= $n90_td){
		$n90 = $_;
		last;
	}
}

$n50 = sprintf ("%.0f", $n50);
$n50 = commify($n50);

$n75 = sprintf ("%.0f", $n75);
$n75 = commify($n75);

$n90 = sprintf ("%.0f", $n90);
$n90 = commify($n90);

# print metrics
$sum_len = commify($sum_len);
print SUM "### Metrics:\n";
print SUM "Number of contigs:\t".$number_of_seq."\n";
print SUM "Total length:\t".$sum_len."\n";
print SUM "Longest contig:\t".$longest."\n";
print SUM "Shortest contig:\t".$shortest."\n";
print SUM "Median length:\t".$median."\n";
print SUM "Average length:\t".$average_len."\n";
print SUM "N50 = $n50\n";
print SUM "N75 = $n75\n";
print SUM "N90 = $n90\n";

###### Checking sequences for telomeres ########################################

my %telo_metrics;

# revc pattern
$pattern = uc($pattern);
my $revc = reverse($pattern);
$revc =~ tr/ATGCatgc/TACGtacg/;

# iterate through sequences
foreach my $seq (sort(keys %sequences)){

	my $start = substr($sequences{$seq}, 0, $window);
	my $end = substr($sequences{$seq}, -$window, $window);

	printseq($seq, $start, $end, \*CHK);

}

##### Priting telomere data to summary #########################################

# Header for telomere data
print SUM "\n\n### Contig\tLength";
print SUM "\t# of pattern (start)\t# of rev pattern (start)";
print SUM "\t# of pattern (end)\t# of rev pattern (end)";
print SUM "\tStatus";
print SUM "\n";

foreach my $contig (sort (keys %telo_metrics)){

	print SUM $contig;
	my $clength = length($sequences{$contig});
	$clength = commify($clength);
	print SUM "\t".$clength;

	my $cf_start = $telo_metrics{$contig}{'f_pattern_start'};
	my $cr_start = $telo_metrics{$contig}{'r_pattern_start'};
	my $cf_end = $telo_metrics{$contig}{'f_pattern_end'};
	my $cr_end = $telo_metrics{$contig}{'r_pattern_end'};

	print SUM "\t".$cf_start."\t".$cr_start;
	print SUM "\t".$cf_end."\t".$cr_end;

	my $telosum = 0;
	if (($cf_start >= $min) or ($cr_start >= $min)){
		$telosum++;
	}
	if (($cf_end >= $min) or ($cr_end >= $min)){
		$telosum++;
	}

	if ($telosum == 2){
		print SUM "\t".'complete chromosome';
	}
	elsif ($telosum == 1){
		print SUM "\t".'partial chromosome with 1 telomere end';
	}
	else {
		print SUM "\t".'contig with no chromosome end attached';
	}

	print SUM "\n";

}

##### Subroutine(s) ############################################################

sub printseq {

	my $contig = $_[0];
	my $start = $_[1];
	my $end = $_[2];
	my $FH1 = $_[3];

	my $count = 0;
	# my @sequences = ($start,$end);
	for my $seq ($start,$end){

		$count++;
		if ($count == 1){
			print CHK "### $contig - START; ";
		}
		else {
			print CHK "### $contig - END; ";
		}

		## Count the number of telomere repeats (forward + reverse)
		my $tmp_f = $seq;
		my @matches_f = ($tmp_f =~ /$pattern/g);
		my $count_f = scalar(@matches_f);

		my $tmp_r = $seq;
		my @matches_r = ($tmp_r =~ /$revc/g);
		my $count_r = scalar(@matches_r);

		print $FH1 "$pattern: $count_f; $revc: $count_r\n";

		if ($count == 1){
			$telo_metrics{$contig}{'f_pattern_start'} = $count_f;
			$telo_metrics{$contig}{'r_pattern_start'} = $count_r;
		}
		else {
			$telo_metrics{$contig}{'f_pattern_end'} = $count_f;
			$telo_metrics{$contig}{'r_pattern_end'} = $count_r;
		}

		## Print substrings (useful to visualize patterns)
		my @seq = unpack("(A60)*", $seq);
		foreach (@seq){
			print $FH1 "$_\n";
		}
		print $FH1 "\n";

	}

}

sub commify {
	my $text = reverse $_[0];
	$text =~ s/(\d{3})(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text;
}