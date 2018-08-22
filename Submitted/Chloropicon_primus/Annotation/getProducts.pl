#!/usr/bin/perl
## Generates a product list from the NCBI faa files

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.faa$//;
	open OUT, ">$file.products";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+)\s+(.*)\s+\[/){
			my $locus = $1;
			my $product = $2;
			print OUT "$locus\t$product\n";
		}
	}
}