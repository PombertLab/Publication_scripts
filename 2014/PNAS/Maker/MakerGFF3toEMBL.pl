#!/usr/bin/perl

## Pombert Lab, IIT, 2014

## Converts the MAKER GFF3 output to preliminary EMBL files for inspection. Writes predicted proteins to a separate FASTA file with the .prot extension.
## NOTE: Requires the  GFF3 outputs (*.gff) and the corresponding fasta inputs (*.fsa) in the same folder

## BUG: Proteins output is buggy, Augustus and Maker work differently. Need to fix this.
## BUG: Locus_tags are out of order. This is because Maker works that way (hashes). Better to implement this separately.

use strict;
use warnings;
use Bio::SeqIO;

my $usage = 'USAGE = GFF3toEMBL *.gff';

die $usage unless @ARGV;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.gff$//;
	open DNA, "<$file.fsa";
	open OUT, ">$file.embl";
	open PROT, ">$file.prot";
	
	### Creating a single DNA string for protein translation
	my $DNAseq= undef;
	while (my $dna = <DNA>){
		chomp $dna;
		if ($dna =~ /^>/){next;}
		else{$DNAseq.= $dna;}
	}
	my $DNAsequence = lc($DNAseq); ## Changing to lower case to fit with the translation hash
	my $contig_length = length($DNAsequence); ## Calculating the contig size
	
	### Init hashes, arrays and values
	my %gene = (); ## Initialize gene hash
	my %exon = (); ## Initialize exon hash
	my %mRNA = (); ## Initialize mRNA hash ### STILL NEEDS TO BE IMPLEMENTED 
	my %intron = (); ## Initialize intron hash
	my %strands = (); ## Initialize strandedness hash
	my @todo = (); ## Initialize array of things to do
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
	my $geneid = 0; ## Initialize gene id
	my $CDS_counter = 0;
	my $contig = undef;
	### Filling the hashes and arrays ###
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^contig-(\d+)\tmaker\tgene\t(\d+)\t(\d+)\t.\t([+-])/){
			$contig =$1;
			my $start = $2;
			my $end = $3;
			my $strandedness = $4;
			$geneid += 10;
			$CDS_counter = 0;
			push (@{$gene{$geneid}}, $start, $end);
			push (@todo, $geneid);
			$strands{$geneid} = $strandedness;
		}
		#~ elsif ($line =~ /\S+\tAUGUSTUS\tintron\t(\d+)\t(\d+)/){
			#~ my $start = $1;
			#~ my $end = $2;
			#~ my $strandedness = $3;
			#~ push (@{$intron{$geneid}}, $1, $2);
		#~ }
		elsif ($line =~ /^contig-\d+\tmaker\tCDS\t(\d+)\t(\d+)\t.\t([+-])\t(\d)/){
			my $start = $1;
			my $end = $2;
			my $strandedness = $3;
			if (($strandedness eq '+') && ( $CDS_counter == 0)){
				push (@{$exon{$geneid}}, $start, $end);
				$CDS_counter++;
			}
			else{
				push (@{$exon{$geneid}}, $start, $end);
			}
		}
	}
	
	### Working on the gene list. Using the @todo array so that we don't have to sort the hash ouput ###
	while (my $list = shift@todo){
		if ($strands{$list} eq '+'){ ## Looking for positive strandedness
			
			### Work on gene features
			print OUT 'FT   gene             '."$gene{$list}[0]".'..'."$gene{$list}[1]"."\n";
			print OUT 'FT                   /locus_tag=DI09_'."$contig".'g'."$list"."\n";
				
			### Work on CDS features
			my $cds_count = scalar(@{$exon{$list}});
			my $end = ($cds_count - 1);
			my $num = ($cds_count - 2);
			my $stopcodon = $exon{$list}[$end];
					
			if (scalar(@{$exon{$list}}) == 2){ ## Verify if we have a single exon
				print OUT 'FT   CDS             '."$exon{$list}[0]..$stopcodon\n";
				print PROT ">$file".'_'."$list\n";
				my $start = (($exon{$list}[0]) - 1);
				my $stop = (($exon{$list}[1]) - 4);
				my $protein = undef;
				for(my $i = $start; $i < ($stop - 2); $i += 3){
					my $codon = substr($DNAsequence, $i, 3);
					$protein .=  $aa{$codon};
				}
				print PROT "$protein\n";
			}
			elsif (scalar(@{$exon{$list}}) > 2){ ## Verify if we have more than one exon
				print OUT 'FT   CDS             '."join($exon{$list}[0]..";
				print PROT ">$file".'_'."$list\n";
				foreach my $count (1..$num){
					if ($count % 2){print OUT "$exon{$list}[$count],";} # Working on odd numbers
					else{print OUT "$exon{$list}[$count]..";} # Working on even numbers
				}
				print OUT "$stopcodon)\n";
	
				my $mRNA = undef;
				my $protein = undef;
				my $tmp1 = undef;
				my $tmp2 = undef;
				foreach my $subs (0..$end){
					if ($subs % 2){# Working on odd numbers
						$tmp2 = ($exon{$list}[$subs]-1);
						$mRNA .= substr($DNAsequence, $tmp1, (($tmp2-$tmp1)+1));
					}
					else{	# Working on even numbers
						$tmp1 = ($exon{$list}[$subs]-1);
					}
				}
				for(my $i = 0; $i < (length($mRNA) - 2); $i += 3){
					my $codon = substr($mRNA, $i, 3);
					$protein .=  $aa{$codon};				
				}
				print PROT "$protein\n";
	
				### Work on intron features
				#~ if (scalar @{$intron{$list}} > 0){
					#~ my $introncount = scalar(@{$intron{$list}});
					#~ my $num = ($introncount - 1);
					#~ foreach my $count (0..$num){
						#~ if ($count % 2){# Working on odd numbers
							#~ print OUT "$intron{$list}[$count]\n";
						#~ }
						#~ else{	# Working on even numbers
							#~ print OUT 'FT   intron          '."$intron{$list}[$count]..";
						#~ }
					#~ }
				#~ }
			}

			### Work on mRNA features ### STILL NEEDS TO BE IMPLEMENTED
		}
		elsif ($strands{$list} eq '-'){ ## Looking for negative strandedness
			
			### Work on gene features
			print OUT 'FT   gene             complement('."$gene{$list}[0]".'..'."$gene{$list}[1]".')'."\n";
			print OUT 'FT                   /locus_tag=DI09_'."$contig".'g'."$list"."\n";
			
			### Work on CDS features
			my $cds_count = scalar(@{$exon{$list}});
			my $end = ($cds_count - 1);
			my $num = ($cds_count - 2);
			my $secs = ($cds_count - 3); ## value for second exon
			my $stopcodon =  $exon{$list}[0]; 
			
			if (scalar(@{$exon{$list}}) == 2){ ## Verify if we have a single exon
				my $start = ($exon{$list}[1] ); 
				my $stop = (($exon{$list}[0]) - 1);
				
				print OUT 'FT   CDS             complement('."$stopcodon..$start)\n";
				print PROT ">$file".'_'."$list\n";
								
				my $protein = undef;
				my $revmRNA =  substr($DNAsequence, $stop, ($start-$stop));
				my $mRNA = reverse($revmRNA);
				$mRNA =~ tr/ATGCRYSWKMBDHVatgcryswkmbdhv/TACGYRWSMKVHDBtacgyrwsmkvhdb/;
				for(my $i = 0; $i < (length($mRNA) - 2); $i += 3){
					my $codon = substr($mRNA, $i, 3);
					$protein .=  $aa{$codon};				
				}
				print PROT "$protein\n";
			}
			elsif (scalar(@{$exon{$list}}) > 2){
				
				my $start = $exon{$list}[$end];
				
				print OUT 'FT   CDS             complement('."join($stopcodon..";
				print PROT ">$file".'_'."$list\n";
				
				foreach my $count (1..$num){
					if ($count % 2){print OUT "$exon{$list}[$count],";} # Working on odd numbers
					else{print OUT "$exon{$list}[$count]..";}	# Working on even numbers
				}
				print OUT "$start))\n";
				
				my $revmRNA = undef;
				my $protein = undef;
				my $tmp1 = undef;
				my $tmp2 = undef;
				foreach my $subs (0..$secs){
					if ($subs % 2){# Working on odd numbers
						$tmp2 = ($exon{$list}[$subs]-1);
						$revmRNA .= substr($DNAsequence, $tmp1, (($tmp2-$tmp1)+1));
					}
					else{	# Working on even numbers
						$tmp1 = ($exon{$list}[$subs]-1);
					}
				}
				foreach my $subset ($num..$end){
					if ($subset % 2){# Working on odd numbers
						$tmp2 = ($exon{$list}[$subset]-1);
						$revmRNA .= substr($DNAsequence, $tmp1, (($tmp2-$tmp1)+1));
					}
					else{	# Working on even numbers
						$tmp1 = ($exon{$list}[$subset]-1);
					}
				}
				my $mRNA = reverse($revmRNA);
				$mRNA =~ tr/ATGCRYSWKMBDHVatgcryswkmbdhv/TACGYRWSMKVHDBtacgyrwsmkvhdb/;
				for(my $i = 0; $i < (length($mRNA) - 2); $i += 3){
					my $codon = substr($mRNA, $i, 3);
					$protein .=  $aa{$codon};				
				}
				print PROT "$protein\n";
				
				### Work on intron features
				#~ if (scalar @{$intron{$list}} > 0){
					#~ my $introncount = scalar(@{$intron{$list}});
					#~ my $num = ($introncount - 1);
					#~ foreach my $count (0..$num){
						#~ if ($count % 2){print OUT "$intron{$list}[$count])\n";} # Working on odd numbers
						#~ else{print OUT 'FT   intron          complement('."$intron{$list}[$count]..";} # Working on even numbers
					#~ }
				#~ }
			}
		}
	}
	### Converting Fasta input to EMBL with BioPerl SeqIO
	my $in = Bio::SeqIO->new(-file => "$file.fsa", -format => "fasta");
	my $out = Bio::SeqIO->new(-file => ">$file.tmp1", -format => "embl");
	while (my $seq = $in->next_seq()){
	$out->write_seq($seq);
	}
	### Adding formatted sequence to the EMBL file
	open SEQ, "<$file.tmp1";
		while (my $read = <SEQ>){
			if ($read =~ /^AC|ID|XX|DE|FH/){next;}
			else{print OUT "$read"};
		}
	close SEQ;
	system "rm $file.tmp1";
}
close IN;
close OUT;
