#!/usr/bin/perl

## Pombert Lab, IIT, 2014
## Version 1.2: Added a verification check for codon_start features
## Writes predicted proteins and mRNAs to separate FASTA files with the .prot and .mRNA extensions.
## NOTE: Requires the EMBL and the corresponding fasta inputs (*.fsa) in the same folder

use strict;
use warnings;
use Bio::SeqIO;

my $usage = 'USAGE = EMBLtoPROT *.embl';

die $usage unless @ARGV;

sub numSort {
		if ($a < $b) { return -1; }
		elsif ($a == $b) { return 0;}
		elsif ($a > $b) { return 1; }
}

my $locus_tag = undef;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.embl$//;
	open DNA, "<$file.fsa";
	open PROT, ">$file.prot";
	open MRNA, ">$file.mRNA";
	
	### Creating a single DNA string for protein translation
	my $DNAseq= undef;
	while (my $dna = <DNA>){
		chomp $dna;
		if ($dna =~ /^>/){next;}
		else{$DNAseq.= $dna;}
	}
	my $DNAsequence = lc($DNAseq); ## Changing to lower case to fit with the translation hash
	my $contig_length = length($DNAsequence); ## Calculating the contig size
	
	my %aa = ('tca'=>'S','tcc'=>'S','tcg'=>'S','tct'=>'S','tcy'=>'S','tcr'=>'S','tcw'=>'S','tcs'=>'S','tck'=>'S','tcm'=>'S','tcb'=>'S','tcd'=>'S','tch'=>'S','tcv'=>'S','tcn'=>'S',
	'ttc'=>'F','ttt'=>'F','tty'=>'F','tta'=>'L','ttg'=>'L','ttr'=>'L','tac'=>'Y','tat'=>'Y','tay'=>'Y','taa'=>'X','tag'=>'X','tga'=>'X','tgc'=>'C','tgt'=>'C','tgy'=>'C',
	'tgg'=>'W','cta'=>'L','ctc'=>'L','ctg'=>'L','ctt'=>'L','cty'=>'L','ctr'=>'L','cts'=>'L','ctw'=>'L','ctk'=>'L','ctm'=>'L','ctb'=>'L','ctd'=>'L','cth'=>'L','ctv'=>'L','ctn'=>'L',
	'cca'=>'P','ccc'=>'P','ccg'=>'P','cct'=>'P','ccr'=>'P','ccy'=>'P','ccs'=>'P','ccw'=>'P','cck'=>'P','ccm'=>'P','ccb'=>'P','ccd'=>'P','cch'=>'P','ccv'=>'P','ccn'=>'P',
	'cac'=>'H','cat'=>'H','cay'=>'H','caa'=>'Q','cag'=>'Q','car'=>'Q','ata'=>'I','atc'=>'I','att'=>'I','atw'=>'I','aty'=>'I','ath'=>'I','atg'=>'M','aac'=>'N','aat'=>'N','aay'=>'N',
	'cga'=>'R','cgc'=>'R','cgg'=>'R','cgt'=>'R','cgr'=>'R','cgy'=>'R','cgs'=>'R','cgw'=>'R','cgk'=>'R','cgm'=>'R','cgb'=>'R','cgd'=>'R','cgh'=>'R','cgv'=>'R','cgn'=>'R',
	'aca'=>'T','acc'=>'T','acg'=>'T','act'=>'T','acr'=>'T','acy'=>'T','acs'=>'T','acw'=>'T','ack'=>'T','acm'=>'T','acb'=>'T','acd'=>'T','ach'=>'T','acv'=>'T','acn'=>'T',
	'aaa'=>'K','aag'=>'K','aar'=>'K','agc'=>'S','agt'=>'S','agy'=>'S','aga'=>'R','agg'=>'R','agr'=>'R','gac'=>'D','gat'=>'D','gay'=>'D','gaa'=>'E','gag'=>'E','gar'=>'E',
	'gta'=>'V','gtc'=>'V','gtg'=>'V','gtt'=>'V','gtr'=>'V','gty'=>'V','gts'=>'V','gtw'=>'V','gtk'=>'V','gtm'=>'V','gtb'=>'V','gtd'=>'V','gth'=>'V','gtv'=>'V','gtn'=>'V',
	'gca'=>'A','gcc'=>'A','gcg'=>'A','gct'=>'A','gcr'=>'A','gcy'=>'A','gcs'=>'A','gcw'=>'A','gck'=>'A','gcm'=>'A','gcb'=>'A','gcd'=>'A','gch'=>'A','gcv'=>'A','gcn'=>'A',
	'gga'=>'G','ggc'=>'G','ggg'=>'G','ggt'=>'G','ggr'=>'G','ggy'=>'G','ggs'=>'G','ggw'=>'G','ggk'=>'G','ggm'=>'G','ggb'=>'G','ggd'=>'G','ggh'=>'G','ggv'=>'G','ggn'=>'G',
	);
	
	### Create hashes to adjust for codon_start features 
	my %coord = ();
	my %codon_start = ();
	my @feat = ();
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /FT                   \/locus_tag="(\DI09_\d+p\d+)"/){$locus_tag = $1; push(@feat, $locus_tag);}
		elsif ($line =~ /FT   CDS/){$coord{$locus_tag} = $line;}
		elsif ($line =~ /FT\s+\/codon_start=(\d+)/){$codon_start{$locus_tag} = $1;}
	}
	
	my $protein = undef;
	my $mRNA = undef;
	while (my $locus = shift@feat){
		chomp $locus;
		my $line = $coord{$locus}; ## Fetch corresponding coordinates line
		### Adjusting codon_start feature, if required
		my $cdn_start = undef;
		if (exists $codon_start{$locus}){$cdn_start = $codon_start{$locus};}
		else{$cdn_start = 1;}
		$cdn_start--;
		### End of adjustment	
		if ($line =~ /FT   CDS             (\d+)..(\d+)/){ ## Forward, single exon
			my $start = $1;
			my $stop = $2;
			$start+=$cdn_start;
			my $lg = ($stop - $start +1);
			$start--;
			$stop--;
			$mRNA = substr($DNAsequence, $start, $lg);
			for(my $i = $start; $i < ($stop - 2); $i += 3){
				my $codon = substr($DNAsequence, $i, 3);
				if (exists  $aa{$codon}){$protein .=  $aa{$codon};}
				else {$protein .= 'X';}
			}
			print MRNA ">$locus\n";
			print MRNA "$mRNA\n";
			print PROT ">$locus\n";
			print PROT "$protein\n";
			$protein=undef;
			$mRNA=undef;
		}
		elsif ($line =~ /FT   CDS             join\((.*)\)/){ ## Forward, multiple exon
			my @array = split(',',$1);
			my @start = ();
			my @stop = ();
			while (my $segment = shift@array){
				chomp $segment;
				if ($segment =~ /(\d+)..(\d+)/){
					my $strt = $1;
					my $stp = $2;
					$strt--;
					$stp--;
					push (@start, $strt);
					push (@stop, $stp);
				}
			}
			$start[0]+=$cdn_start;
			my $asize = scalar(@start);
			my $num = ($asize -1);
			foreach my $subs (0..$num){
					$mRNA .= substr($DNAsequence, $start[$subs], (($stop[$subs]-$start[$subs])+1));
			}
			for(my $i = 0; $i < (length($mRNA) - 5); $i += 3){
				my $codon = substr($mRNA, $i, 3);
				if (exists  $aa{$codon}){$protein .=  $aa{$codon};}
				else {$protein .= 'X';}			
			}
			print MRNA ">$locus\n";
			print MRNA "$mRNA\n";
			print PROT ">$locus\n";
			print PROT "$protein\n";
			$protein=undef;
			$mRNA=undef;
		}
		elsif ($line =~ /FT   CDS             complement\((\d+)..(\d+)\)/){ ## Reverse, single exon
			my $start = $2;
			my $stop = $1;
			$start-=$cdn_start;
			$start--;
			$stop--;
			my $revmRNA =  substr($DNAsequence, ($stop), ($start-$stop+1));
			$mRNA = reverse($revmRNA);
			$mRNA =~ tr/ATGCRYSWKMBDHVatgcryswkmbdhv/TACGYRWSMKVHDBtacgyrwsmkvhdb/;
			for(my $i = 0; $i < (length($mRNA) - 5); $i += 3){
				my $codon = substr($mRNA, $i, 3);
				if (exists  $aa{$codon}){$protein .=  $aa{$codon};}
				else {$protein .= 'X';}
			}
			print MRNA ">$locus\n";
			print MRNA "$mRNA\n";
			print PROT ">$locus\n";
			print PROT "$protein\n";
			$protein=undef;
			$mRNA=undef;
		}
		elsif ($line =~ /FT   CDS             complement\(join\((.*)/){ ## Reverse, mutiple exon
			my @array = split(',',$1);
			my @start = ();
			my @stop = ();
			while (my $segment = shift@array){
				chomp $segment;
				if ($segment =~ /(\d+)..(\d+)/){
					my $strt = $1;
					my $stp = $2;
					$strt--;
					$stp--;
					unshift (@start, $strt);
					unshift (@stop, $stp);
				}
			}
			
			my @sstart = sort numSort @start;
			my @sstop = sort numSort @stop;
			my $asize = scalar(@start);
			my $num = ($asize -1);
			$sstop[$num]-=$cdn_start;
			my $revmRNA = undef;
			foreach my $subs (0..$num){
					$revmRNA .= substr($DNAsequence, $sstart[$subs], (($sstop[$subs]-$sstart[$subs])+1));
			}
			$mRNA = reverse($revmRNA);
			$mRNA =~ tr/ATGCRYSWKMBDHVatgcryswkmbdhv/TACGYRWSMKVHDBtacgyrwsmkvhdb/;
			for(my $i = 0; $i < (length($mRNA) - 5); $i += 3){
				my $codon = substr($mRNA, $i, 3);
				if (exists  $aa{$codon}){$protein .=  $aa{$codon};}
				else {$protein .= 'X';}			
			}
			print MRNA ">$locus\n";
			print MRNA "$mRNA\n";
			print PROT ">$locus\n";
			print PROT "$protein\n";
			$protein=undef;
			$mRNA = undef;
		}
	}
}
close IN;
close PROT;
