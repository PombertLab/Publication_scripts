#!/usr/bin/env perl
# Pombert Lab, 2023

my $name = 'gbff_parser.pl';
my $version = '0.1a';
my $updated = '2023-10-09';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);
use PerlIO::gzip;

my $usage = <<"USAGE";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Parses GenBank GBFF files for CDS, rRNA and tRNA features, then
            summarizes keys metrics

COMMAND     ${name} -i *.gbff -o ./

OPTIONS:
-i (--input)    GenBank GBFF input file(s) # GZIP supported; file type determined by file extension
-o (--outdir)   Output directory [Default: ./]
USAGE

die "\n$usage\n" unless @ARGV;

my @input_files;
my $outdir = './';
GetOptions(
    'i|input=s@{1,}' => \@input_files,
    'o|outdir=s' => \$outdir
);

unless (-d $outdir){
    make_path($outdir,{mode=>0755}) or die("Can't create output directory $outdir: $!\n");
}

foreach my $input_file (@input_files){

    my $file_name = basename($input_file);
    my @file_data = split('\.',$file_name);
    my $file_prefix = $file_data[0];
    my $diamond = "<";
    my $ext = $file_data[-1];

    if ($file_data[-1] eq 'gz'){ 
        $diamond = "<:gzip";
        $ext = $file_data[-2];
    }

    print "Parsing $file_name\n";

    my $comprehensive = $outdir.'/'.$file_prefix.'.details.txt';
    my $summary = $outdir.'/'.$file_prefix.'.summary.txt';

    open OUT, '>', $comprehensive or die "Can't write to output file: $comprehensive\n";
    open MET, '>', $summary or die "Can't write to output file: $summary\n";
    open GBK, $diamond, "$input_file" or die "Can't open input file $input_file: $!\n";

    ### Database per file
    my %features;
    my %metrics;
    my %sequences;

    my $feature;
     my $sum_contig_length;

    my $locus;
    my $product;

    my $incomplete_product_name;
    my $incomplete_sections;

    my $contig;
    my $strand;
    my $origin;

    while (my $line = <GBK>){

        chomp $line;

        if ($line =~ /^LOCUS\s+(\S+)\s+(\d+) bp/ ){
            $contig = $1;
            $metrics{'contig_num'}++;
            $sum_contig_length += $2;
        }

        ## Grabbing sequences
        if ($line =~ /^ORIGIN/){
            $origin = 1;
        }
        elsif ($line =~ /^\/\//){
            $origin = undef;
            next;
        }
        elsif (($origin) and ($line !~ /^ORIGIN/)){
            my $seq = $line;
            $seq =~ s/\d//g;
            $seq =~ s/\s//g;
            $seq = uc($seq);
            $sequences{$locus} .= $seq;
        }

        ## Locus tags
        if ($line =~ /^\s{21}\/locus_tag="(.*)"/){
            $locus = $1;
            $locus =~ s/\W/\_/g;
        }

        # Grabbing CDS, rRNA and tRNA features: using mRNA 
        # instead of CDS because some GBFF were lacking 
        # the CDS feature despite having the associated mRNA tag

        if ($line =~ /^\s{5}(\w+)/){
            $feature = $1;
        }

        if ($line =~ /^\s{5}(mRNA|rRNA|tRNA)\s+(?:complement)*\(*<*(\d+)\.\.>*(\d+)/){

            $features{$locus}{'contig'} = $contig;

            $features{$locus}{'feature'} = $1;
            push (@{$features{$locus}{'exons'}}, "$2\.\.$3");

            if ($line =~ /complement/){
                $strand = '-';
            }
            else {
                $strand = '+';
            }

            $features{$locus}{'strand'} = $strand;
            $strand = undef;

        }

        elsif ($line =~ /^\s{5}(mRNA|rRNA|tRNA)\s+(?:complement\()*(?:order\()*(?:join\()*(.*)\)$/){

            $features{$locus}{'contig'} = $contig;

            $features{$locus}{'feature'} = $1;
            my $section = $2;
            $section =~ s/complement//g;
            $section =~ s/join//g;
            $section =~ s/\(//g;
            $section =~ s/\)//g;
            $section =~ s/\<//g;
            $section =~ s/\>//g;

            my @sections = split(",",$section);
            push (@{$features{$locus}{'exons'}}, @sections);

            if ($line =~ /complement/){
                $strand = '-';
            }
            else {
                $strand = '+';
            }

            $features{$locus}{'strand'} = $strand;
            $strand = undef;

        }

        elsif ($line =~ /^\s{5}(mRNA|rRNA|tRNA)\s+(?:complement\()*(?:order\()*(?:join\()*(.*)[^)]$/){

            $features{$locus}{'feature'} = $1;
            $features{$locus}{'contig'} = $contig;
            $incomplete_sections = $2;

            if ($line =~ /complement/){
                $strand = '-';
            }
            else {
                $strand = '+';
            }

            $features{$locus}{'strand'} = $strand;
            $strand = undef;

        }
        elsif (defined $incomplete_sections){

            if ($line =~ /^\s{21}(.*)\)$/){

                $incomplete_sections .= ','.$1;

                $incomplete_sections =~ s/complement//g;
                $incomplete_sections =~ s/join//g;
                $incomplete_sections =~ s/\(//g;
                $incomplete_sections =~ s/\)//g;
                $incomplete_sections =~ s/\>//g;
                $incomplete_sections =~ s/\<//g;
                $incomplete_sections =~ s/,,/,/g;

                my @sections = split(",", $incomplete_sections);
                push (@{$features{$locus}{'exons'}}, @sections);

                $incomplete_sections = undef;

            }
            elsif ($line =~ /^\s{21}(.*)/){
                $incomplete_sections .= ','.$1;
            }

        }

        ## Checking for product names, including accross multiple lines
        if (defined $feature){

            if ($feature eq ('mRNA') or $feature eq ('rRNA') or $feature eq ('tRNA')){

                if ($line =~ /^\s{21}\/product="(.*)"/){
                    $product = $1;
                    $features{$locus}{'product'} = $product;
                }

                elsif ($line =~ /^\s{21}\/product="(.*)/){
                    $incomplete_product_name = $1;
                }

                elsif ((defined $incomplete_product_name) && ($line =~ /^\s{21}(.*)\"$/)){

                    $incomplete_product_name .= ' ';
                    $incomplete_product_name .= $1;
                    $product = $incomplete_product_name;
                    $incomplete_product_name = undef;

                    $features{$locus}{'product'} = $product;
                    $product = undef;

                }

                elsif ((defined $incomplete_product_name) && ($line =~ /^\s{21}(.*)/)){
                    $incomplete_product_name .= ' ';
                    $incomplete_product_name .= $1;
                }

            }

        }

    }


    ### Calculating metrics
    print OUT '# Contig'."\t";
    print OUT 'Locus tag'."\t";
    print OUT 'Type'."\t";
    print OUT 'Strand'."\t";
    print OUT 'Start'."\t";
    print OUT 'End'."\t";
    print OUT 'Exons'."\t";
    print OUT 'Ex. length'."\t";
    print OUT 'Introns'."\t";
    print OUT 'In. length'."\t";
    print OUT 'Product'."\n";

    foreach my $key (sort (keys %features)){

        $metrics{'total'}++;

        my $feature = $features{$key}{'feature'};
        if ($feature eq 'mRNA'){
            $feature = 'CDS';
        }

        my $strand = $features{$key}{'strand'};
        $metrics{$feature}++;
        $metrics{$strand}++;

        ## Locus information
        print OUT $features{$key}{'contig'}."\t";
        print OUT $key."\t";
        print OUT $feature."\t";
        print OUT $strand."\t";

        ## Exons
        my @exons = @{$features{$key}{'exons'}};
        my $num_exons = scalar(@exons);
        $metrics{'exon_total'} += $num_exons;

        my $first_ex = $exons[0];
        my $last_ex = $exons[-1];

        my @ex_first = split('\.\.', $first_ex);
        my @ex_last = split('\.\.', $last_ex);

        my $ex_start = $ex_first[0];
        my $ex_end = $ex_last[-1];

        print OUT $ex_start."\t";
        print OUT $ex_end."\t";
        print OUT $num_exons."\t";

        my $exon_length = 0;
        for my $exon (@exons){
            my ($start, $end) = split('\.\.', $exon);
            my $length = $end - $start + 1;
            $exon_length += $length;
            $metrics{'exon_total_length'} += $length;
        }

        print OUT $exon_length."\t";

        ## Introns (if any)
        my $num_introns = $num_exons - 1;
        $metrics{'intron_total'} += $num_introns;
        print OUT $num_introns."\t";

        my $intron_length = 0;
        if ($num_introns >= 1){

            $metrics{'genes_with_introns'}++;

            for my $x (1..$num_introns){

                my $slice_1 = $exons[$x-1];
                my $slice_2 = $exons[$x];

                my @ex_1 = split('\.\.', $slice_1);
                my @ex_2 = split('\.\.', $slice_2);

                my $start = $ex_1[-1];
                my $end = $ex_2[0];
                my $i_length = $end - $start + 1;

                $intron_length += $i_length;
                $metrics{'intron_total_length'} += $i_length;

            }

        }

        print OUT $intron_length."\t";
        print OUT $features{$key}{'product'}."\n";

    }

    ## Writing FASTA out
    my $fasta_out = $outdir.'/'.$file_prefix.'.fasta';
    open FASTA, '>', $fasta_out or die "Can't write to output file: $fasta_out\n";
    foreach my $contig (sort (keys %sequences)){
        print FASTA '>'.$contig."\n";
        my @seq = unpack ("(A60)*", $sequences{$contig});
		while (my $seq = shift@seq){
			print FASTA "$seq\n";
		}
    }
    close FASTA;

    ## Quick metrics summary
    my ($data_entry) = fileparse($input_file);

    print MET '# '.$data_entry."\n";

    print MET 'Contigs:'."\t".$metrics{'contig_num'}."\n";
    print MET 'Length (bp):'."\t".$sum_contig_length."\n";
    print MET 'Features:'."\t".$metrics{'total'}."\n";
    print MET '# of CDS:'."\t".$metrics{'CDS'}."\n";

    if (!defined $metrics{'rRNA'}){
        $metrics{'rRNA'} = 0;
    }
    if (!defined $metrics{'tRNA'}){
        $metrics{'tRNA'} = 0;
    }
    print MET '# of rRNA:'."\t".$metrics{'rRNA'}."\n";
    print MET '# of tRNA:'."\t".$metrics{'tRNA'}."\n";

    my $average_gene_density = ($metrics{'total'} * 1000) / $sum_contig_length;
    $average_gene_density = sprintf("%.3f", $average_gene_density);
    print MET 'Av. gene/kb:'."\t".$average_gene_density."\n";

    my $average_exon_per_gene = $metrics{'exon_total'}/$metrics{'total'};
    $average_exon_per_gene = sprintf("%.2f", $average_exon_per_gene);
    print MET 'Av. exon/gene:'."\t".$average_exon_per_gene."\n";

    my $average_intron_per_gene = $metrics{'intron_total'}/$metrics{'total'};
    $average_intron_per_gene = sprintf("%.2f", $average_intron_per_gene);
    print MET 'Av. intron/gene:'."\t".$average_intron_per_gene."\n";

    my $average_exon_length = $metrics{'exon_total_length'}/$metrics{'exon_total'};
    $average_exon_length = sprintf("%.0f", $average_exon_length);
    print MET 'Av. exon length:'."\t".$average_exon_length."\n";

    my $average_intron_length = $metrics{'intron_total_length'}/$metrics{'intron_total'};
    $average_intron_length = sprintf("%.0f", $average_intron_length);
    print MET 'Av. intron length:'."\t".$average_intron_length."\n";

    my $percent_genes_introns = ($metrics{'genes_with_introns'}/$metrics{'total'}) * 100;
    $percent_genes_introns = sprintf("%.2f", $percent_genes_introns);
    print MET '# of genes with introns:'."\t".$metrics{'genes_with_introns'}."\n";
    print MET '% of genes with introns:'."\t".$percent_genes_introns.'%'."\n";

}

