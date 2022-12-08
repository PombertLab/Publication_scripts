#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2019
my $version = '0.6';
my $name = 'run_pilon.pl';
my $updated = '2021-04-07';

use strict; use warnings; use Getopt::Long qw(GetOptions); use File::Basename;

## Defining options
my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Runs pilon read correction in paired-end mode for X iterations (stops automatically 
		if Pilon no longer makes changes to the consensus)

REQUIREMENTS	SAMtools; minimap2; pilon

EXAMPLE		run_pilon.pl \\
		  -f polished_genome_flye.fasta \\
		  -pe1 EEint_S1_L001_R1_001.fastq \\
		  -pe2 EEint_S1_L001_R2_001.fastq \\
		  -r 10 \\
		  -p /opt/pilon/pilon-1.22.jar \\
		  -t 64

OPTIONS:
-f (--fasta)	Genome assembly in FASTA format
-pe1		Illumina R1 file in FASTQ format
-pe2		Illumina R2 file in FASTQ format
-r (--rounds)	Number of pilon iterations [Default: 10]
-p (--pilon)	Path to pilon jar file [Default: /opt/pilon-1.24/pilon-1.24.jar]
-t (--threads)	Number of threads [Default: 16]
-g (--Gb)	Amount of RAM (in Gb) to allocate to pilon jar file [Default: 16]
-o (--outdir)	Output directory [Default: ./]
USAGE
die "\n$usage\n" unless @ARGV;
my @commands = @ARGV;

my $fasta;
my $pe1;
my $pe2;
my $rounds = 10;
my $pilon = '/opt/pilon-1.24/pilon-1.24.jar';
my $threads = 16;
my $ram = 16;
my $outdir = './';
GetOptions(
	'f|fasta=s' => \$fasta,
	'pe1=s' => \$pe1,
	'pe2=s' => \$pe2,
	'r|rounds=i' => \$rounds,
	'p|pilon=s' => \$pilon,
	't|threads=i' => \$threads,
	'g|Gb=i' => \$ram,
	'o|outdir=s' => \$outdir
);

## Checking output directory and log file
unless (-d $outdir) {
	mkdir ($outdir,0755) or die "Can't create $outdir: $!\n";
}
open LOG, ">", "$outdir/pilon.log" or die "Can't create $outdir/pilon.log: $!\n";
my $timestamp = `date`;
chomp $timestamp;
print LOG "Ran pilon on $timestamp with $name version $version\n\n";
print LOG "COMMAND: $name @commands\n";

## Basename + variables
my ($file, $path) = fileparse($fasta);
$file =~ s/\.\w+$//;

my $sam_file = "$outdir/round1.sam";
my $bam_unsorted = "$outdir/round1.bam";
my $bam_sorted = "$outdir/round1.sorted.bam";
my $round = "$outdir/round1";

### Creation bam files and run them through pilon
# BAM round no. 1
system "minimap2 -t $threads -R \@RG\\\\tID:$pe1\\\\tSM:$fasta -ax sr $fasta $pe1 $pe2 1> $sam_file";
system "samtools view -@ $threads -bS $sam_file -o $bam_unsorted";
system "samtools sort -@ $threads -o $bam_sorted  $bam_unsorted";
system "rm $bam_unsorted $sam_file"; ## Discarding unsorted BAM/SAM files
system "samtools index $bam_sorted ";

# Pilon round no. 1
system "java \\
  -Xmx${ram}G \\
  -jar $pilon \\
  --genome $fasta \\
  --frags $bam_sorted \\
  --output round1 \\
  --outdir $round \\
  --changes";

# Cleanup
system "sed 's/_pilon//g' $round/round1.fasta > $round/round1.clean.fasta";
system "mv $round/round1.clean.fasta $round/round1.fasta";

# Subsequent rounds
my $prev;
my $num;
if ($rounds >= 2){
	for $num (2..($rounds)){
		$prev = $num - 1;

		## Checking if changes were made in the previous iteration. If not, stop.
		my $filename = "${outdir}/round${prev}/round${prev}.changes";
		my $size = -s $filename;
		last if ($size == 0);
		system "mv ${outdir}/round${prev}.sorted.bam ${outdir}/round${prev}.sorted.bam.bai $outdir/round${prev}/";

		## Naming files again
		$sam_file = "$outdir/round${num}.sam";
		$bam_unsorted = "$outdir/round${num}.bam";
		$bam_sorted = "$outdir/round${num}.sorted.bam";
		$round = "$outdir/round${num}";

		### Running additional pilon iterations if changes were made in the previous iteration
		# BAM round 2+
		system "minimap2 -t $threads -R \@RG\\\\tID:$pe1\\\\tSM:$fasta -ax sr ${outdir}/round${prev}/round${prev}.fasta $pe1 $pe2 1> $sam_file";
		system "samtools view -@ $threads -bS $sam_file -o $bam_unsorted";
		system "samtools sort -@ $threads -o $bam_sorted $sam_file";
		system "rm $bam_unsorted $sam_file"; ## Discarding unsorted BAM/SAM files
		system "samtools index $bam_sorted";

		# Pilon round 2+
		system "java \\
		  -Xmx${ram}G \\
		  -jar $pilon \\
		  --genome ${outdir}/round${prev}/round${prev}.fasta \\
		  --frags $bam_sorted \\
		  --output round${num} \\
		  --outdir $round \\
		  --changes";

		# Cleanup
		system "sed 's/_pilon//g' $round/round${num}.fasta > $round/round${num}.clean.fasta";
		system "mv $round/round${num}.clean.fasta $round/round${num}.fasta";
	}
}
print "No changes after Pilon iteration # $prev. Stopping subsequent iterations.\n";

# Final cleanup
system "mv ${outdir}/round${prev}.sorted.bam ${outdir}/round${prev}.sorted.bam.bai $outdir/round${prev}/";
system "cp $outdir/round${prev}/round${prev}.fasta $outdir/$file.pilon.fasta";
