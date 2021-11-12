#!/usr/bin/perl
## Pombert lab, 2021
my $name = 'parallel_COFACTOR.pl';
my $version = '0.1b';
my $updated = '2021-05-25';

use strict; use warnings; use Getopt::Long qw(GetOptions);
use File::Basename; use threads; use threads::shared;
use Cwd 'abs_path';

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Runs parralel instances of COFACTOR

NOTE	If path too long, COFACTOR chokes...

COMMAND		${name} \\
		  -c /media/Data_3/I-TASSER/I-TASSER5.1/I-TASSERmod/runCOFACTOR.pl \\
		  -l /media/Data_3/I-TASSER/I-TASSER5.1/Libraries/ \\
		  -p /media/Data_3/jpombert/Microsporidia/PDBs/ \\
		  -o /media/Data_3/jpombert/Microsporidia/COFACTOR_RESULTS

OPTIONS:
-c (--cof)		Location of runCOFACTOR.pl
-l (--lib)		Locations of I-TASSER libraries
-p (--pdb)		Location of PDB files to query
-o (--out)		Output directory [Default: ./]
-t (--threads)	Number of cpu threads to use [Default: 10] ## i.e. runs n processes in parallel
OPTIONS
die "\n$usage\n" unless @ARGV;

my $cofactor;
my $libdir;
my $pdbs;
my $outdir = './';
my $threads = 10;
GetOptions(
	'c|cof=s' => \$cofactor,
	'l|lib=s' => \$libdir,
	'p|pdb=s' => \$pdbs,
	'o|out=s' => \$outdir,
	't|threads=i' => \$threads
);

## Add PDBs from directory to array
opendir (my $dh, $pdbs) or die "Can't opendir $pdbs: $!\n";
my @pdbs;
while (readdir $dh) {
	my $fname = $_;
	if ($fname =~ /.pdb$/){
		push (@pdbs, $fname);
		print "$fname\n";
	}
}
closedir $dh;

## Working on output directory
unless (-d $outdir){
	mkdir ($outdir,0755) or die "Can't create $outdir: $!\n";
	$outdir = abs_path($outdir);
}
else {
	$outdir = abs_path($outdir);
}

## multi-threading
my @threads = initThreads(); ## Initialize # of threads specified
my @files :shared = @pdbs; ## Copying the array into a shared list for multithreading (use threads::shared;)
@files = sort (@files);
for(@threads){ $_ = threads->create(\&exe); }	# Tell threads run the exe sub
for(@threads){ $_ -> join(); }	# Run until threads are done
exit;

## Subroutines
sub initThreads{ 
	# An array to place our threads in
	my @initThreads;
	for (my $i = 1; $i <= $threads; $i++){ push(@initThreads,$i); }
	return @initThreads;
}

sub exe{

	my $id = threads->tid();  #Get the thread id. Allows each thread to be identified.

	## Iterating on pdb files
	while (my $pdb = shift @files){
		
		print "\nThread $id working on $pdb...\n";
		my ($pdb_name) = $pdb =~ /^(\S+)\.pdb/;
				
		## Copying + changing to subdir to prevent overwriting temporary files used by COFACTOR
		my $tmp_dir = "$outdir/$pdb_name";
		unless (-d $tmp_dir){ 
			mkdir ($tmp_dir,0755) or die "Can't create $tmp_dir: $!\n";
		}
		system "cp $pdbs/$pdb_name.pdb $tmp_dir/";

		chdir "$tmp_dir";

		# Running COFACTOR
		system "$cofactor \\
			-model $pdb \\
			-protname $pdb_name \\
			-libdir $libdir \\
			-datadir $tmp_dir \\
			-LBS false \\
			-EC false
		";

	}
	threads->exit();
}
