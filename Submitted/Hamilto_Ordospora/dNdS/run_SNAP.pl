#!/usr/bin/perl
## Pombert Lab, 2018
## Runs SNAP on aligned sequences in fasta format
## v1.0; now implements multithreading

use strict; use warnings;
use threads; use threads::shared;
use Getopt::Long qw(GetOptions);

my $usage = <<'OPTIONS';

USAGE = run_SNAP.pl -t 10 -s SNAP_mod.pl -a *.macse
OPTIONS:
-s	## Modified SNAP script
-a	## Alignment in multifasta format
-t	## Threads [Default: 10]

OPTIONS
die $usage unless @ARGV;

my $threads = 10;
my @fasta;
my $snap;
GetOptions(
	't=i' => \$threads,
	'a=s@{1,}' => \@fasta,
	's=s' => \$snap
);
my @files :shared = @fasta;
my @threads = initThreads();
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
	while (my $fasta = shift@files){
		print "Thread $id working on $fasta...\n";
		open IN, "<$fasta";
		$fasta =~ s/.fasta$//; $fasta =~ s/.fa$//; $fasta =~ s/.fsa$//; $fasta =~ s/.macse$//;
		open OUT, ">$fasta.sequence";
		my %sequences = ();
		my $seq; my @seqs = ();
		while (my $line = <IN>){
			chomp $line;
			if ($line =~ /^>(\S+)/){$seq = $1; push(@seqs, $seq);}
			else{
				$line = uc $line;
				$line =~ s/!/-/g;	# replaces bang (!) for gaps (-) so that SNAP can work on the output files
				$sequences{$seq} .= $line;
			}
		}
		while (my $tmp = shift @seqs){print "found $tmp\n";print OUT "$tmp\t$sequences{$tmp}\n";}
		system "perl $snap $fasta.sequence";
	}	
	threads->exit();
}