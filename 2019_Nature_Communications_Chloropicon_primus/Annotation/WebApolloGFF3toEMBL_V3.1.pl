#!/usr/bin/perl

## Pombert Lab, IIT, 2016

## Converts WebApollo GFF3 files to EMBL files and writes annotated proteins and mRNAs to a separate FASTA files with the .prot and .mRNA extensions.
## Generates locus tags automatically based on the provided prefix.
## NOTE: Requires the GFF3 outputs (*.gff3) and the corresponding fasta inputs (*.fsa) in the same folder
## NOTE: Requires the BioPerl module Bio::SeqIO

use strict;
use warnings;
use Bio::SeqIO;

my $usage = 'USAGE = WebApolloGFF3toEMBL.pl LOCUS_TAG_PREFIX *.gff3';
die $usage unless @ARGV;

my $locus_tag_prefix = shift@ARGV;
my $locus_id = 10;
my $contig_number = 0;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.gff3$//;
	open DNA, "<$file.fsa";
	open OUT, ">$file.embl";
	open PROT, ">$file.prot";
	open MRNA, ">$file.mRNA";
	$contig_number++;
	
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
	my $geneid = 0; ## Initialize gene id
	my $CDS_counter = 0; ## Initialize CDS counter
	my %gene = (); ## Initialize gene hash
	my %exon = (); ## Initialize exon hash
	my %strands = (); ## Initialize strandedness hash
	my @todo = (); ## Initialize array of things to do
	my %aa = ('tca'=>'S','tcc'=>'S','tcg'=>'S','tct'=>'S','tcy'=>'S','tcr'=>'S','tcw'=>'S','tcs'=>'S','tck'=>'S','tcm'=>'S','tcb'=>'S','tcd'=>'S','tch'=>'S','tcv'=>'S','tcn'=>'S',
	'ttc'=>'F','ttt'=>'F','tty'=>'F','tta'=>'L','ttg'=>'L','ttr'=>'L','tac'=>'Y','tat'=>'Y','tay'=>'Y','taa'=>'*','tag'=>'*','tga'=>'*','tgc'=>'C','tgt'=>'C','tgy'=>'C',
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
		
	### Filling the hashes and arrays ###
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /\S+\t\S+\tgene\t(\d+)\t(\d+)\t\S+\t([-+])/){
			my $start = $1;
			my $end = $2;
			my $strandedness = $3;
			$geneid = $start;
			push (@{$gene{$geneid}}, $start, $end);
			push (@todo, $geneid);
			$strands{$geneid} = $strandedness;
			$CDS_counter = 0;
		}
		elsif ($line =~ /\S+\t\S+\tCDS\t(\d+)\t(\d+)\t\S+\t([-+])\t(\d)/){
			my $start = $1;
			my $end = $2;
			my $strandedness = $3;
			my $phase = $4;
			if ($strandedness eq '+'){
				push (@{$exon{$geneid}}, $start, $end);
				$CDS_counter++;
			}
			elsif ($strandedness eq '-'){
				push (@{$exon{$geneid}}, $end, $start);
				$CDS_counter++;
			}
		}
	}
	
	### Working on the gene list. Using the @todo array so that we don't have to sort the hash ouput ###
	my @sorted_todo = sort {$a <=> $b} @todo;
	while (my $list = shift@sorted_todo){
		if ($strands{$list} eq '+'){ ## Looking for positive strandedness
			
			### Work on gene features
			my $locus_number = sprintf("%05d", $locus_id);
			my $contig = sprintf("%02d", $contig_number);
			print OUT 'FT   gene             '."$gene{$list}[0]".'..'."$gene{$list}[1]"."\n";
			print OUT 'FT                   /locus_tag="'."$locus_tag_prefix".'_'."$contig".'p'."$locus_number".'"'."\n";
			$locus_id += 10;
				
			### Work on CDS features
			my $cds_count = scalar(@{$exon{$list}});
			my $end = ($cds_count - 1);
			my $num = ($cds_count - 2);
			my $stopcodon = $exon{$list}[$end];			
			
			if (scalar(@{$exon{$list}}) == 2){ ## Verify if we have a single exon
				print OUT 'FT   CDS             '."$exon{$list}[0]..$stopcodon\n";
				print OUT 'FT                   /locus_tag="'."$locus_tag_prefix".'_'."$contig".'p'."$locus_number".'"'."\n";
				print PROT ">$locus_tag_prefix".'_'."$contig".'p'."$locus_number\n";
				print MRNA ">$locus_tag_prefix".'_'."$contig".'p'."$locus_number\n";
				my $start = (($exon{$list}[0]));
				my $stop = (($exon{$list}[1]));
				my $protein = undef;
				my $mRNA =  substr($DNAsequence, $start-1, (($stop-$start)+1));
				for(my $i = 0; $i < (length($mRNA) - 5); $i += 3){ ## -2 if stop codon (*) desired, -5 if not
					my $codon = substr($mRNA, $i, 3);
					$protein .=  $aa{$codon};				
				}
				my @PROTEIN = unpack ("(A60)*", $protein);
				while (my $tmp = shift@PROTEIN){
					print PROT "$tmp\n";
				}
				my @RNA = unpack ("(A60)*", $mRNA);
				while (my $tmp = shift@RNA){
					print MRNA "$tmp\n";
				}
			}
			elsif (scalar(@{$exon{$list}}) > 2){ ## Verify if we have more than one exon
				print OUT 'FT   CDS             '."join($exon{$list}[0]..";
				print PROT ">$locus_tag_prefix".'_'."$contig".'p'."$locus_number\n";
				print MRNA ">$locus_tag_prefix".'_'."$contig".'p'."$locus_number\n";
				foreach my $count (1..$num){
					if ($count % 2){print OUT "$exon{$list}[$count],";} # Working on odd numbers
					else{print OUT "$exon{$list}[$count]..";} # Working on even numbers
				}
				print OUT "$stopcodon)\n";
				print OUT 'FT                   /locus_tag="'."$locus_tag_prefix".'_'."$contig".'p'."$locus_number".'"'."\n";
				
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
				for(my $i = 0; $i < (length($mRNA) - 5); $i += 3){ ## -2 if stop codon (*) desired, -5 if not
					my $codon = substr($mRNA, $i, 3);
					$protein .=  $aa{$codon};				
				}
				my @PROTEIN = unpack ("(A60)*", $protein);
				while (my $tmp = shift@PROTEIN){
					print PROT "$tmp\n";
				}
				my @RNA = unpack ("(A60)*", $mRNA);
				while (my $tmp = shift@RNA){
					print MRNA "$tmp\n";
				}
			}
		}
		if ($strands{$list} eq '-'){ ## Looking for positive strandedness
			
			### Work on gene features
			my $locus_number = sprintf("%05d", $locus_id);
			my $contig = sprintf("%02d", $contig_number);
			print OUT 'FT   gene             complement('."$gene{$list}[0]".'..'."$gene{$list}[1]".')'."\n";
			print OUT 'FT                   /locus_tag="'."$locus_tag_prefix".'_'."$contig".'p'."$locus_number".'"'."\n";
			$locus_id += 10;
				
			### Work on CDS features
			my $cds_count = scalar(@{$exon{$list}});
			my $end = ($cds_count - 1);
			my $num = ($cds_count - 2);
			my @reversed = reverse(@{$exon{$list}}); ## Reversing the list of exons
			my $stopcodon = $reversed[$end];	
			
			if (scalar(@reversed) == 2){ ## Verify if we have a single exon
				print OUT 'FT   CDS             complement('."$reversed[0]..$reversed[1]".')'."\n";
				print OUT 'FT                   /locus_tag="'."$locus_tag_prefix".'_'."$contig".'p'."$locus_number".'"'."\n";
				print PROT ">$locus_tag_prefix".'_'."$contig".'p'."$locus_number\n";
				print MRNA ">$locus_tag_prefix".'_'."$contig".'p'."$locus_number\n";
				my $start = (($reversed[0]));
				my $stop = (($reversed[1]));
				my $protein = undef;
				my $revmRNA =  substr($DNAsequence, $start-1, (($stop-$start)+1));
				my $mRNA = reverse($revmRNA);
				$mRNA =~ tr/ATGCRYSWKMBDHVatgcryswkmbdhv/TACGYRWSMKVHDBtacgyrwsmkvhdb/;
				for(my $i = 0; $i < (length($mRNA) - 5); $i += 3){ ## -2 if stop codon (*) desired, -5 if not
					my $codon = substr($mRNA, $i, 3);
					$protein .=  $aa{$codon};				
				}
				my @PROTEIN = unpack ("(A60)*", $protein);
				while (my $tmp = shift@PROTEIN){
					print PROT "$tmp\n";
				}
				my @RNA = unpack ("(A60)*", $mRNA);
				while (my $tmp = shift@RNA){
					print MRNA "$tmp\n";
				}				
			}
			elsif (scalar(@reversed) > 2){ ## Verify if we have more than one exon
				print OUT 'FT   CDS             complement('."join($reversed[0]..";
				print PROT ">$locus_tag_prefix".'_'."$contig".'p'."$locus_number\n";
				print MRNA ">$locus_tag_prefix".'_'."$contig".'p'."$locus_number\n";
				foreach my $count (1..$num){
					if ($count % 2){print OUT "$reversed[$count],";} # Working on odd numbers
					else{print OUT "$reversed[$count]..";} # Working on even numbers
				}
				print OUT "$stopcodon)\n";
				print OUT 'FT                   /locus_tag="'."$locus_tag_prefix".'_'."$contig".'p'."$locus_number".'"'."\n";
	
				my $revmRNA = undef;
				my $protein = undef;
				my $tmp1 = undef;
				my $tmp2 = undef;
				foreach my $subs (0..$end){
					if ($subs % 2){# Working on odd numbers
						$tmp2 = ($reversed[$subs]-1);
						$revmRNA .= substr($DNAsequence, $tmp1, (($tmp2-$tmp1)+1));
					}
					else{	# Working on even numbers
						$tmp1 = ($reversed[$subs]-1);
					}
				}
				my $mRNA = reverse($revmRNA);
				$mRNA =~ tr/ATGCRYSWKMBDHVatgcryswkmbdhv/TACGYRWSMKVHDBtacgyrwsmkvhdb/;
				for(my $i = 0; $i < (length($mRNA) - 5); $i += 3){ ## -2 if stop codon (*) desired, -5 if not
					my $codon = substr($mRNA, $i, 3);
					$protein .=  $aa{$codon};				
				}
				my @PROTEIN = unpack ("(A60)*", $protein);
				while (my $tmp = shift@PROTEIN){
					print PROT "$tmp\n";
				}
				my @RNA = unpack ("(A60)*", $mRNA);
				while (my $tmp = shift@RNA){
					print MRNA "$tmp\n";
				}
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
