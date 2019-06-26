#!/usr/bin/perl
## Pombert Lab, 2018
my $name = 'run_macse.pl';
my $version = '1.0';

use strict; use warnings; use Getopt::Long qw(GetOptions);
use threads; use threads::shared;

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Runs MACSE alignments on multifasta files 
USAGE		run_macse.pl -t 10 -f *.fasta -v /opt/MACSE/macse_v2.01.jar -p alignSequences -g 1

OPTIONS:
-t	## Threads [Default: 10]
-f	## Files in multifasta format
-v	## Path to MACSE jar file [Default: /opt/MACSE/macse_v2.01.jar]
-p	## MACSE program [Default: alignSequences]
-g	## NCBI genetic code [Default: 1]

OPTIONS
die $usage unless @ARGV;

my $threads = 10;
my @fasta;
my $var = '/opt/MACSE/macse_v2.01.jar';
my $program = 'alignSequences';
my $gc = '1';
GetOptions(
	't=i' => \$threads,
	'f=s@{1,}' => \@fasta,
	'v=s' => \$var,
	'p=s' => \$program,
	'g=i' => \$gc
);

my @threads = initThreads();
my @files :shared = @fasta;	## Copying the array into a shared list for multithreading (use threads::shared;)
for(@threads){$_ = threads->create(\&exe);}	# Tell threads run the exe sub
for(@threads){$_->join();}	# Run until threads are done
exit;

## subroutines
sub initThreads{ # An array to place our threads in
	my @initThreads;
	for(my $i = 1;$i<=$threads;$i++){push(@initThreads,$i);}
	return @initThreads;
}
sub exe{
	my $id = threads->tid();  #Get the thread id. Allows each thread to be identified.
	while (my $fasta = shift @files) {
		$fasta =~ s/.fasta//;
		print "Thread $id aligning $fasta...\n";
		system "java -jar $var -prog $program -seq $fasta.fasta -out_NT $fasta.NT.macse -out_AA $fasta.AA.macse -gc_def $gc";
	}	
	threads->exit();
}
