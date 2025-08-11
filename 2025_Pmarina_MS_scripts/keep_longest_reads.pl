#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2020
my $version = '0.9b';
my $name = 'keep_longest_reads.pl';
my $updated = '2022-10-21';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use PerlIO::gzip;
use File::Basename;
use File::Path qw(make_path);

################################################################################
# Usage definition
################################################################################

my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	KLR (keep longest reads) calculates metrics for FASTQ file(s)
		and/or parses them to keep the longest reads either by minimum size
		or by desired sequencing depth (useful for large Nanopore or PacBio
		datasets).

COMMAND LINE EXAMPLES
Minimum read length:		${name} -i *.fastq -o ./KLR -m 10000
Desired sequencing depth:	${name} -i file.fastq -o ./KLR -d 100 -s 3000000
Metrics only:			${name} -i *.fastq -o ./KLR -m 10000 -x -t -j -h 'Long read data'

I/O OPTIONS:
-i (--input)	Input file(s) in FASTQ format
-o (--outdir)	Output directory [Default: ./]
-p (--prefix)	Desired prefix for log, tsv and json files [Default: klr_metrics]
-t (--tsv)	Metrics summary in TSV format
-j (--json)	Create .json file for multiQC
-h (--head)	JSON section name for multiQC [Default: Long read data]
-q (--qcid)	ID name for multiQC JSON file [Default: KLR]

PARSING OPTIONS:
-x (--metrics)	Calculate metrics only, do not create read subset
-m (--minimum)	Minimum read length to keep
-d (--depth)	Desired sequencing depth (requires estimated genome size: -s)
-s (--size)	Expected genome size

NOTES:
- The -m and -d options are mutually exclusive
- The -d option should only be used with single FASTQ files or from genomes with
  very similar expected sizes
USAGE

die "\n$usage\n" unless@ARGV;

################################################################################
# Command line options
################################################################################

my @commands = @ARGV;
my @fastq;
my $outdir = './';
my $prefix = 'klr_metrics';
my $tsv;
my $json;
my $json_header = 'Long read data';
my $qcid = 'KLR';
my $metrics;
my $min;
my $depth;
my $genome_size;

GetOptions(
	# i/o
	'i|input=s@{1,}' => \@fastq,
	'o|outdir=s' => \$outdir,
	'p|prefix=s' => \$prefix,
	't|tsv' => \$tsv,
	'j|json' => \$json,
	'h|head=s' => \$json_header,
	'q|qcid=s' => \$qcid,
	# parsing
	'x|metrics'	=> \$metrics,
	'm|minimum=i' => \$min,
	'd|depth=i' => \$depth,
	's|size=i' => \$genome_size
);

################################################################################
# Outdir + logfile
################################################################################

unless (-d $outdir) {
	make_path( $outdir, { mode => 0755 } )  or die "Can't create $outdir: $!\n";
}

my $stime = `date`;
chomp $stime;

my $logfile = $prefix.'.log';
open LOG, ">", "$outdir/$logfile" or die "Can't create $outdir/$logfile: $!\n";
print LOG "COMMAND: $name @commands\n";
print LOG "Started on $stime\n";

################################################################################
# Working on FASTQ file(s)
################################################################################

my %metrics_data;
my $basename;

while (my $fastq = shift@fastq){

	print "\nWorking on $fastq...\n";

	my $gzip = '';
	if ($fastq =~ /.gz$/){ $gzip = ':gzip'; }
	open FASTQ1, "<$gzip", "$fastq" or die "Can't open $fastq: $!\n";
	open FASTQ2, "<$gzip", "$fastq" or die "Can't open $fastq: $!\n";

	## Grabbing basename from files (minus extensions...)
	$basename = fileparse($fastq);
	$basename =~ s/\.gz$//;
	$basename =~ s/\.fastq$//;
	$basename =~ s/\.fq$//;

	## parsing by minimum length
	if ($min or $metrics){

		## Output file
		my $filename;
		if ($min){
			my $div = ($min/1000);
			$filename = $outdir."/".$basename.'_'.$div."k.fastq";
		}
		unless ($metrics){
			open OUT, ">", "$filename" or die "Can't create $filename: $!\n";
		}

		## making sure that $minimum_len is at least 1 to prevent div by 0
		my $minimum_len;
		if ($min){ $minimum_len = $min; }
		else {
			$minimum_len = 1;
		}

		my @lengths;
		my @subset;
		my %reads;
		my $count = 0;
		my $read;

		## Working on FASTQ file
		while (my $line = <FASTQ1>){
			chomp $line;

			## Read name
			if ($count == 0){
				$read = $line;
				$reads{$read}[0] = $line;
				$count++;
			}

			## Read sequence
			elsif ($count == 1){ 
				$reads{$read}[1] = $line;
				my $len = length($line);
				push (@lengths, $len);
				$count++;
			}

			## Skipping + sign
			elsif ($count == 2){ $count++; }

			## Quality score 
			elsif ($count == 3){
				$reads{$read}[2] = $line;
				if (length($reads{$read}[1]) >= $minimum_len){
					my $keep = length($reads{$read}[1]);
					push (@subset, $keep);
					unless ($metrics){
						print OUT "$reads{$read}[0]\n";
						print OUT "$reads{$read}[1]\n";
						print OUT '+'."\n";
						print OUT "$reads{$read}[2]\n";
					}
				}
				## Clearing read db to minimize memory usage
				$count = 0; %reads = (); 
			}
		}

		if ($fastq =~ /.gz$/){ binmode FASTQ1, ":gzip(none)"; }	

		## Calculating metrics
		metrics($fastq, 'full', \@lengths); ## Full
		if ($min){
			if ($min != 1){ ## min == 1 is the same as the full dataset...
				metrics($filename, $min, \@subset); ## subset
			}
		}
		unless ($metrics){
			close OUT;
		}

	}

	elsif ($depth){

		## Output file
		my $filename = $outdir."/".$basename.'_'.$min."x.fastq";
		unless ($metrics){
			open OUT, ">", "$filename" or die "Can't create $filename: $!\n";
		}

		my @lengths;
		my $count = 0;
		## Doing two independent passes to avoid loading the full FASTQ file in memory
		## Pass 1 - Calculating metrics
		while (my $line = <FASTQ1>){
			chomp $line;
			if ($count == 0){ $count++; }
			elsif ($count == 1){
				my $len = length($line);
				push (@lengths, $len);
				$count++;
			}
			elsif ($count == 2){ $count++; }
			elsif ($count == 3){ $count=0; }
		}

		if ($fastq =~ /.gz$/){ binmode FASTQ1, ":gzip(none)"; }	

		@lengths = sort {$b <=> $a} @lengths; ## sort by size, from largest to smallest
		my $len_threshold;
		my $sum;
		foreach (@lengths){
			if ($sum <= ($depth*$genome_size)){ $len_threshold = $_; }
			$sum += $_;
		}
		$len_threshold = $len_threshold;
		print "\nMinimum read size for ${depth}X sequencing depth at estimated genome size of $genome_size bp = $len_threshold bp\n";
		print "Saving reads of at least $len_threshold bp to $filename\n\n";

		## Pass 2 - Parsing reads
		$count = 0;
		my $read;
		my %reads;
		my @subset;

		while (my $line = <FASTQ2>){
			chomp $line;

			## Read name
			if ($count == 0){
				$read = $line;
				$reads{$read}[0] = $line;
				$count++;
			}

			## Read sequence
			elsif ($count == 1){ 
				$reads{$read}[1] = $line;
				$count++;
			}

			## Skipping + sign
			elsif ($count == 2){ $count++; }

			## Quality score
			elsif ($count == 3){
				$reads{$read}[2] = $line; 
				if (length($reads{$read}[1]) >= $len_threshold){
					my $keep = length($reads{$read}[1]);
					push (@subset, $keep);
					print OUT "$reads{$read}[0]\n";
					print OUT "$reads{$read}[1]\n";
					print OUT '+'."\n";
					print OUT "$reads{$read}[2]\n";
				}
				## Clearing read db to minimize memory usage
				$count = 0;
				%reads = ();
			}
		}

		if ($fastq =~ /.gz$/){ binmode FASTQ2, ":gzip(none)"; }

		## Calculating metrics
		metrics($fastq, 'full', \@lengths); ## Full
		metrics($filename, $depth, \@subset); ## Subset
		close OUT;
	}

	close FASTQ1;
	close FASTQ2;

}

################################################################################
# TSV metrics table
################################################################################

if ($tsv){
	my $tsv_file = $outdir.'/'.$prefix.'.tsv';
	open TSV, ">", "$tsv_file" or die "Can't create $tsv_file: $!\n";

	# Header
	foreach my $key (sort (keys %metrics_data)){
		print TSV "\t$key";
	}
	print TSV "\n";

	# Metrics
	tsv(\*TSV, 'Number of reads', 'nreads');
	tsv(\*TSV, 'Total number of bases', 'total');
	tsv(\*TSV, 'Longest read', 'longest');
	tsv(\*TSV, 'Shortest read', 'shortest');
	tsv(\*TSV, 'Average read size', 'average');
	tsv(\*TSV, 'Median read size', 'median');
	tsv(\*TSV, 'N50', 'n50');
	tsv(\*TSV, 'N75', 'n75');
	tsv(\*TSV, 'N90', 'n90');

}

################################################################################
# .json metrics file for MultiQC
################################################################################

if ($json){

	my $json_file = $outdir.'/'.$prefix.'_mqc.json';
	open MQC, ">", "$json_file" or die "Can't create $json_file: $!\n";

	my $header = <<"	HEAD";
	"id": "$qcid",
	"section_name": "$json_header",
	"description": "Metrics calculated with keep_longest_reads.pl (https://github.com/PombertLab).",
	"plot_type": "table",
	"pconfig": {
		"id": "$qcid",
		"format": "{:,.0f}"
	},
	HEAD

	print MQC "\{\n";
	print MQC "$header";
	print MQC "\t\"data\": {\n";
	my $comma_counter = 0;

	foreach my $key (sort (keys %metrics_data) ){
		print MQC "\t\t\"$key\": {\n";
		print MQC "\t\t\t\"Number of reads\": $metrics_data{$key}{'nreads'},\n";
		print MQC "\t\t\t\"Total number of bases\": $metrics_data{$key}{'total'},\n";
		print MQC "\t\t\t\"Longest read\": $metrics_data{$key}{'longest'},\n";
		print MQC "\t\t\t\"Shortest read\": $metrics_data{$key}{'shortest'},\n";
		print MQC "\t\t\t\"Average read size\": $metrics_data{$key}{'average'},\n";
		print MQC "\t\t\t\"Median read size\": $metrics_data{$key}{'median'},\n";
		print MQC "\t\t\t\"N50\": $metrics_data{$key}{'n50'},\n";
		print MQC "\t\t\t\"N75\": $metrics_data{$key}{'n75'},\n";
		print MQC "\t\t\t\"N90\": $metrics_data{$key}{'n90'}\n";
		unless ($comma_counter == scalar (keys %metrics_data) - 1){
			print MQC "\t\t},\n";
		}
		else {
			print MQC "\t\t}\n";
		}
		$comma_counter++;
	}

	print MQC "\t}\n}";

}

################################################################################
# Final timestamp
################################################################################

my $etime = `date`;
chomp $etime;

unless ($metrics){
	print LOG "Ended on: $etime\n";
}

################################################################################
# Subroutine(s)
################################################################################

# Prints desired metric in tab-delimited columns (TSV) to filehandle
sub tsv {

	my $fh = $_[0];
	my $column = $_[1];
	my $metric = $_[2];
	print $fh "$column";
	foreach my $key (sort (keys %metrics_data)){
		my $value = $metrics_data{$key}{$metric};
		# Comment out the line below if you want a TSV without commas in numbers
		$value = commify($value);
		print $fh "\t$value";
	}
	print $fh "\n";

}

# Calculates metrics from FASTQ file
sub metrics {

	my @fh = (*LOG, *STDOUT);
	my $file = $_[0];

	## basename
	my $basename_ext = $_[1];
	unless ($basename_ext eq 'full'){
		if ($min){
			my $div = ($basename_ext/1000);
			$basename = $basename." >= ".$div.'k';
		}
		elsif ($depth){
			$basename = $basename." ".$basename_ext.'x';
		}
		
	}

	# number of reads in dataset
	my @reads = @{$_[2]};
	my $num_reads = scalar @reads;
	my @len = sort {$b <=> $a} (@reads); ## sort by size; from largest to smallest

	$metrics_data{$basename}{'nreads'} = $num_reads;
	my $nreads = commify($num_reads);
	foreach (@fh){
		print $_ "\n## Metrics for dataset $file\n\n";
		print $_ "Number of reads: $nreads\n";
	}

	# median
	my $median;
	my $median_pos = $num_reads/2;
	if ($median_pos =~ /^\d+$/){ $median = $len[$median_pos]; }
	else {
		my $med1 = int($median_pos);
		my $med2 = $med1 + 1;
		$median = (($len[$med1] + $len[$med2])/2);
	}
	$median = sprintf("%.0f", $median);
	$metrics_data{$basename}{'median'} = $median;
	$median = commify($median);

	# sum
	my $sum;
	foreach (@len){ $sum += $_; }
	$metrics_data{$basename}{'total'} = $sum;
	my $fsum = commify($sum);

	# longest
	my $large = sprintf("%.0f", $len[0]);
	$metrics_data{$basename}{'longest'} = $large;
	$large = commify($large);

	# shortest
	my $small = sprintf("%.0f", $len[$#len]);
	$metrics_data{$basename}{'shortest'} = $small;
	$small = commify($small);

	# average
	my $average = sprintf("%.0f", ($sum/$num_reads));
	$metrics_data{$basename}{'average'} = $average;
	$average = commify($average);

	foreach (@fh){
		print $_ "Total number of bases: $fsum\n";
		print $_ "Longest read = $large nt\n";
		print $_ "Shortest read = $small nt\n";
		print $_ "Average read size = $average nt\n";
		print $_ "Median read size = $median nt\n";
	}

	### N50, N75, N90
	# thresholds to reach for N metrics
	my $n50_td = $sum*0.5;
	my $n75_td = $sum*0.75;
	my $n90_td = $sum*0.9;
	# sums to calculate
	my $nsum50 = 0;
	my $nsum75 = 0;
	my $nsum90 = 0;
	# metrics to capture
	my $n50;
	my $n75,
	my $n90;

	foreach (@len){
		$nsum50 += $_;
		if ($nsum50 >= $n50_td){
			$n50 = $_;
			last;
		}
	}
	foreach (@len){
		$nsum75 += $_;
		if ($nsum75 >= $n75_td){
			$n75 = $_;
			last
		}
	}
	foreach (@len){
		$nsum90 += $_;
		if ($nsum90 >= $n90_td){
			$n90 = $_;
			last;
		}
	}

	$n50 = sprintf ("%.0f", $n50);
	$metrics_data{$basename}{'n50'} = $n50;
	$n50 = commify($n50);

	$n75 = sprintf ("%.0f", $n75);
	$metrics_data{$basename}{'n75'} = $n75;
	$n75 = commify($n75);

	$n90 = sprintf ("%.0f", $n90);
	$metrics_data{$basename}{'n90'} = $n90;
	$n90 = commify($n90);

	foreach (@fh){
		print $_ "N50 = $n50 nt\n";
		print $_ "N75 = $n75 nt\n";
		print $_ "N90 = $n90 nt\n";
		print $_ "\n";
	}

}

# Adds commas to numbers; from the Perl Cookbook (O'Reilly)
sub commify {
	my $text = reverse $_[0];
	$text =~ s/(\d{3})(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text;
}