#!/usr/bin/perl

use strict;
use warnings;

## Converts the output of the extractGFF_features.pl script to a TBL format to use with TBL2ASN ##
## NOTE: TBL2ASN doesn't require the various features to be interleaved, so we can output a block of gene features, and then a block of CDS features, and so forth...## 

my $usage = 'perl features_to_TBL.pl *.feat';

die $usage unless @ARGV;

while (my $feat = shift @ARGV) {

	open IN, "<$feat" or die "cannot open $feat";
	$feat =~ s/\.feat$//;
	open OUT, ">$feat.tbl";

	print OUT ">Feature $feat\.fsa\n";
	
	my @genes = ();		## create empty array of genes to be filled later ##
	my @CDS = ();		## create empty array of CDS to be filled later ##
	my @mRNAs = ();		## create empty array of mRNAs to be filled later ##
	my @exons = ();		## create empty array of exons to be filled later ##
	my @introns = ();		## Generate introns later when everything else works ##
	
	
	##############################################################################################################
	
	
	while (my $line = <IN>)  {	## reads the info from the features file and push at the end of the corresponding empty arrays up above ##
		
		chomp $line;
		
		if ($line ~~ /^(\+|\-)\t(match)\t(\d+)\t(\d+)\t(\S+)$/) {
			push (@genes, $line);
		}
		elsif ($line ~~ /^(\+|\-)\t(match_part)\t(\d+)\t(\d+)\t(\S+)$/) {
			push (@CDS, $line);
			push (@exons, $line);
			push (@mRNAs, $line);
		}
	}
	
	
	##############################################################################################################
	
	
	foreach my $genes (@genes) {
		
		if ($genes =~ /^(\+|\-)\t(match)\t(\d+)\t(\d+)\t(\S+)$/) {	## looking for genes entries ##

			my $strandedness = $1;
			my $type = $2;
			my $start = $3;
			my $end = $4;
			my $locus_tag = $5;

				if ($type ~~ 'match' && $strandedness ~~ '+') {	## For genes on the plus strand ##
					print OUT "<$start\t>$end\tgene\n";
					print OUT "\t\t\tlocus_tag\t$locus_tag\n";
				}
				elsif ($type ~~ 'match' && $strandedness ~~ '-') {	## For genes on the minus strand ##
					print OUT "<$end\t>$start\tgene\n";
					print OUT "\t\t\tlocus_tag\t$locus_tag\n";
				}
		}
	}
	
	
	##############################################################################################################
	
	
	my $CDS_tag = undef;
	my $gene_tag = undef;
	
	
	while (my $CDS = shift @CDS) {## Use a until loop with shift @ until @exons is empty ##
			
		if ($CDS ~~ /^(\+)\t(match_part)\t(\d+)\t(\d+)\t(\S+)$/) { ## Let's work on the plus strand first ##
			
			my $start = $3;
			my $end = $4;
			$gene_tag = $5;
		
		
			if (!defined($CDS_tag)) {
				$CDS_tag = $gene_tag;
				print OUT "$start\t$end\tCDS\n";
			}
				
			elsif ($CDS_tag ne $gene_tag) {
				print OUT "\t\t\tlocus_tag\t$CDS_tag\n";
				#~ print OUT "\t\t\tproduct\thypothetical protein\n";		## The tag product is required by TBL2ASN to generate the protein translation in the SQN file ##
				print OUT "\t\t\tprotein_id\tgnl\|ubcbot\|$CDS_tag\n";
				print OUT "\t\t\ttranscript_id\tgnl\|ubcbot\|$CDS_tag\_mRNA\n";
				$CDS_tag = $gene_tag; ## change to the new CDS	
				print OUT "$start\t$end\tCDS\n";		
			}
			
			elsif ($CDS_tag eq $gene_tag) { 						
				print OUT "$start\t$end\n";
			}
				
		}
		
		elsif ($CDS ~~ /^(\-)\t(match_part)\t(\d+)\t(\d+)\t(\S+)$/) { ## Let's work on the minus strand now ##
			
			my $start = $4;
			my $end = $3;
			$gene_tag = $5;
		
		
			if (!defined($CDS_tag)) {
				$CDS_tag = $gene_tag;
				print OUT "$start\t$end\tCDS\n";	
			}
		
			elsif ($CDS_tag ne $gene_tag) {
				print OUT "\t\t\tlocus_tag\t$CDS_tag\n";
				#~ print OUT "\t\t\tproduct\thypothetical protein\n";		## The tag product is required by TBL2ASN to generate the protein translation in the SQN file ##
				print OUT "\t\t\tprotein_id\tgnl\|ubcbot\|$CDS_tag\n";
				print OUT "\t\t\ttranscript_id\tgnl\|ubcbot\|$CDS_tag\_mRNA\n";
				$CDS_tag = $gene_tag;	## change to the new CDS
				print OUT "$start\t$end\tCDS\n";	
			}
			
			elsif ($CDS_tag eq $gene_tag) { 
				print OUT "$start\t$end\n";	
			}
					
		}

	}
	
	if ((scalar(@CDS) == 0) && (scalar(@genes) != 0 )) { ## Print out the remaining info. That should take care of the end of file  ##
		print OUT "\t\t\tlocus_tag\t$CDS_tag\n";
		#~ print OUT "\t\t\tproduct\thypothetical protein\n";		## The tag product is required by TBL2ASN to generate the protein translation in the SQN file ##
		print OUT "\t\t\tprotein_id\tgnl\|ubcbot\|$CDS_tag\n";
		print OUT "\t\t\ttranscript_id\tgnl\|ubcbot\|$CDS_tag\_mRNA\n";
	}
	
	
	##############################################################################################################
	
	
	my $mRNA_tag = undef;
	my $mRNA_gene_tag = undef;
	
	while (my $mRNAs= shift @mRNAs) {## Use a until loop with shift @ until @mRNAs is empty ##
			
		if ($mRNAs ~~ /^(\+)\t(match_part)\t(\d+)\t(\d+)\t(\S+)$/) { ## Let's work on the plus strand first ##
			
			my $start = $3;
			my $end = $4;
			$mRNA_gene_tag = $5;
		
		
			if (!defined($mRNA_tag)) {
				$mRNA_tag = $mRNA_gene_tag;
				print OUT "<$start\t$end\tmRNA\n";
			}
				
			elsif ($mRNA_tag ne $mRNA_gene_tag) {
				print OUT "\t\t\ttranscript_id\tgnl\|ubcbot\|$mRNA_tag\_mRNA\n";
				$mRNA_tag = $mRNA_gene_tag; ## change to the new mRNA	
				print OUT "<$start\t$end\tmRNA\n";		
			}
			
			elsif ($mRNA_tag eq $mRNA_gene_tag) { 						
				print OUT "$start\t$end\n";
			}
				
		}
		
		elsif ($mRNAs ~~ /^(\-)\t(match_part)\t(\d+)\t(\d+)\t(\S+)$/) { ## Let's work on the minus strand now ##
			
			my $start = $4;
			my $end = $3;
			$mRNA_gene_tag = $5;
		
		
			if (!defined($mRNA_tag)) {
				$mRNA_tag = $mRNA_gene_tag;
				print OUT "<$start\t$end\tmRNA\n";	
			}
		
			elsif ($mRNA_tag ne $mRNA_gene_tag) {
				print OUT "\t\t\ttranscript_id\tgnl\|ubcbot\|$mRNA_tag\_mRNA\n";
				$mRNA_tag = $mRNA_gene_tag;	## change to the new mRNA
				print OUT "<$start\t$end\tmRNA\n";	
			}
			
			elsif ($mRNA_tag eq $mRNA_gene_tag) { 
				print OUT "$start\t$end\n";	
			}
					
		}

	}
	
	if (scalar(@mRNAs) == 0 && (scalar(@genes) != 0 )) { ## Print out the remaining info. That should take care of the end of file  ##
		print OUT "\t\t\ttranscript_id\tgnl\|ubcbot\|$mRNA_tag\_mRNA\n";
	}


	##############################################################################################################
	
	
	close IN;
}

exit;
