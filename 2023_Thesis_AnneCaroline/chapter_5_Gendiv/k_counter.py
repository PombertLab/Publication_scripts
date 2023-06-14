#!/usr/bin/env python3

name = 'k_counter.py'
version = '0.4.2'
updated = '2023-04-05'

from sys import argv
from argparse import ArgumentParser
from os import makedirs
from os.path import exists, basename
from re import match
from time import time
from multiprocessing import Pool, cpu_count

usage = f"""
NAME		{name}
VERSION		{version}
UPDATED		{updated}
SYNOPIS		Determines the kmer redundancy of a genome at shifting frames

USAGE		{name} \\
		  -i Eint_2019_11_18_FINAL.fasta \\
		  -o kmer_redundancy \\
		  -t 6 \\
		  -k 5 6 7 8 9 11 13 \\
		  -w 600 \\
		  -s 1

OPTIONS
-i (--input)		Input genome fasta
-o (--outdir)		Output directory name [Default: kmer_redundancy]
-t (--threads)		Number of threads [Default: 1]
-k (--kmers)		Kmer(s) to use [Default: 6]
-m (--mersonly)		Calculate only the kmer redundnacy
-g (--gconly)		Calculate only the GC content
-w (--window)		Window size to check for kmers [Default: 1000]
-s (--slide)		Sliding window size [Default: 1000]
"""

if len(argv) < 2:
	print(usage)
	exit()

GetOptions = ArgumentParser()

GetOptions.add_argument('-i','--input',required=True,nargs="*")
GetOptions.add_argument('-o','--outdir',default="kmer_redundancy")
GetOptions.add_argument('-t','--threads',default=1)
GetOptions.add_argument('-k','--kmers',default=[6],nargs="*")
GetOptions.add_argument('-m','--mersonly',action='store_true',default=False)
GetOptions.add_argument('-g','--gconly',action='store_true',default=False)
GetOptions.add_argument('-w','--window',default=1000)
GetOptions.add_argument('-s','--slide',default=500)

args = GetOptions.parse_args()

in_files = args.input
outdir = args.outdir
basedir = outdir
threads = int(args.threads) if (int(args.threads) <= cpu_count()) else cpu_count()
kmer_sizes = [int(k) for k in args.kmers]
mersonly = args.mersonly
gconly = args.gconly
window_size = int(args.window)
slide = int(args.slide)

status = 'all'
if mersonly and not gconly:
	status = 'mers'
elif gconly and not mersonly:
	status = 'gc'

sequences = {}

for orgs in in_files:
	IN = open(orgs,"r") or exit(f"Cannot read from {orgs}!")
	organism = match("^(\w+)",basename(orgs)).groups()[0]
	basedir = f"{basedir}/{organism}"
	sequences[organism] = {}
	chromosome = ""
	for line in IN:
		line = line.replace("\n","")
		if match("^\>(\S+)",line): 
			chromosome = match("^\>(\S+)",line).groups()[0]
			sequences[organism][chromosome] = {'seq':"", 'kmers':kmer_sizes, 'window_size':window_size,'slide':slide,'chromosome':chromosome,'outdir':f"{outdir}/{organism}/{chromosome}",'status':status}
		else:
			sequences[organism][chromosome]['seq'] += line
	IN.close()

def chromosome_crawl(chromosome_data):

	chromosome = chromosome_data['chromosome']
	outdir = chromosome_data['outdir']
	sequence = chromosome_data['seq']
	length = len(sequence)
	kmer_sizes = chromosome_data['kmers']
	window_size = chromosome_data['window_size']
	slide = chromosome_data['slide']
	status = chromosome_data['status']

	if not exists(outdir):
		makedirs(outdir)

	if status == 'gc' or status == 'all':
		data = ""
		start = 0
		tick = time()
		print(f"Calculating {chromosome} GC content for {int(length-window_size)} windows... This may take a while...")
		while start < length-window_size-1:
			nuc_cont = {"A":0,"C":0,"G":0,"T":0}
			seq = sequence[start:start+window_size]
			for nuc in seq:
				if nuc.upper() in nuc_cont.keys():
					nuc_cont[nuc.upper()] += 1
			data += f"{start}\t{(nuc_cont['C']+nuc_cont['G'])/(sum(nuc_cont.values()))}\n"
			start += slide

		GC_OUT = open(f"{outdir}/{chromosome}.window_{window_size}.slide_{slide}.gc_content","w")
		GC_OUT.write(data)
		GC_OUT.close()
		print(f"Completed GC content calculations for {chromosome}! Duration {time()-tick} seconds!")

	if status == 'mers' or status == 'all':
		for kmer_size in kmer_sizes:
			data = ""
			uk_circos_data = f"#chr START END UNIQUE_KMERS_REDUNDANCY\n"
			k50_circos_data = f"#chr START END UNIQUE_KMERS_REDUNDANCY\n"
			tick = time()
			start = 0
			print(f"Calculating {window_size-kmer_size+1} kmers for {int(length-window_size)} windows on {chromosome}... This may take a while...")

			## Possible kmers per window
			pkpw = window_size-kmer_size+1

			## Possible unique kmers
			puk = 4**kmer_size

			## Expected unique kmer usage
			eku = pkpw/puk

			data += "### WINDOW_LOC_BP\tREDUNDNACY BY UNIQUE_KMERS\tREDUNDANCY BY K50\tREDUNDANCY BY ABOVE AVERAGE KMERS\tAVERAGED REDUNDANCY\n"
			while start < length-window_size:

				kmers = {}
				for begin in range(start,(start+window_size-kmer_size+1)):

					kmer = sequence[begin:begin+kmer_size]
					
					if kmer in kmers.keys():
						kmers[kmer] += 1
					else:
						kmers[kmer] = 1

				## Identified unique kmers found in the window
				iuk = len(kmers.keys())

				## Actual unique kmer usage
				aku = iuk/puk

				# print()
				# print(f"Possible unique kmers (size={kmer_size}): {puk}")
				# print(f"Possible kmers per window (size={kmer_size}): {pkpw}")
				# print(f"Expected usage of unique kmers: {pkpw}/{puk} => {eku}")
				# print(f"Identified unique kmers (size={kmer_size}): {iuk}")
				# print(f"Actual usage of unique kmers: {iuk}/{puk} => {aku}")
				
				## Adjusted possible kmers per window
				apkpw = pkpw

				# If the number of unique kmers found is equal to the potential unique kmers, by default, the calculation
				# will result in a redundancy of 0, even though there may be instances were, given a large sliding window,
				# all possible unique kmers are found, and a set of kmers is expressed in triple the amount of the others.
				# Here we iterate through all possible kmers, reducing a 'layer of redundancy' by subtracting a single
				# instance of a kmer for each possible kmer, and removing a possible kmer if there are no more instances of it.
				# As such, as kmer instances are removed, the possible kmers per window (pkpw) must be adjusted. If at
				# the end, there is only a single kmer that is expressed once more than all the rest, dividing 1 by
				# the pkpw (likely very large) will result in a large redundancy value, even though the kmers were expressed
				# very evenly. The adjusted possible kmers per window (apkpw) is set to intially be equal to the pkpw,
				# and with each layer of redundancy removed, apkpw is reduced by puk, to account for the decrease in the
				# adjusted window size

				if (iuk == puk):
					t_kmers = kmers
					while iuk == puk:
						apkpw -= len(t_kmers.keys())
						remove = []
						for kmer in t_kmers.keys():
							if t_kmers[kmer] == 1:
								remove.append(kmer)
							else:
								t_kmers[kmer] -= 1
						for kmer in remove:
							del t_kmers[kmer]
						iuk = len(t_kmers.keys())
						aku = iuk/puk

				# print(f"Actual usage of unique kmers: {iuk}/{puk} => {aku}")
				
				## Redundancy by unique kmers
				rbuk = (1-(iuk/apkpw)) if apkpw < puk else (1-(iuk/puk))
				# print(f"Redundancy by unique kmers: {rbuk}")

				kmer_values_sorted = sorted(kmers.values(),reverse=True)
				# print(kmer_values_sorted)
				
				threshold = .5*apkpw if apkpw < puk else .5*puk

				## Least number of kmers to produce 50% of kmers in sequence
				k50 = 0
				kmer_total = 0

				while kmer_total < threshold:
					kmer_total += kmer_values_sorted.pop(0)
					k50 += 1

				## Redundancy by K50
				rbk50 =(1-((k50)/(apkpw*.5))) if eku <= 1 else (1-((k50)/(puk*.5)))

				# print(f"Expected K50: {(pkpw*.5) if eku <= 1 else (puk*.5)}")
				# print(f"K50: {k50}")
				# print(f"Redundancy (by K50): {rbk50}")

				## Above average kmers
				aak = sum(1 for i in kmers.values() if i > eku)
				## Redundancy by number of kmers above the expected redundancy of kmers (erok)
				rbaak = (aak/puk) if puk < apkpw else (aak/apkpw)

				# print(f"Number of kmers above expected average: {aak}")
				# print(f"Redundancy (by aak kmers): {rbaak}")

				## Average redundancy
				ar = (rbuk + rbk50 + rbaak)/3

				# print(f"Average: {(rbk50+rbuk+rbaak)/3}")
				# print()

				# exit()

				data += f"{start}\t{rbuk}\t{rbk50}\t{rbaak}\t{ar}\n"
				end = start + window_size - 1
				uk_circos_data += f"{chromosome} {start} {end} {rbuk}\n"
				k50_circos_data += f"{chromosome} {start} {end} {rbk50}\n"
				start += slide
			OUT = open(f"{outdir}/{chromosome}.kmer_{kmer_size}.window_{window_size}.slide_{slide}.kmer_count","w")
			OUT.write(data)
			OUT.close()

			UK_CIRCOS = open(f"{basedir}/kmers.uk.k{kmer_size}.s{slide}.w{window_size}.circos", "a")
			UK_CIRCOS.write(uk_circos_data)
			UK_CIRCOS.close()

			K50_CIRCOS = open(f"{basedir}/kmers.k50.k{kmer_size}.s{slide}.w{window_size}.circos", "a")
			K50_CIRCOS.write(k50_circos_data)
			K50_CIRCOS.close()

			print(f"Completed kmer counting for {chromosome}! Duration {time()-tick} seconds!")


for org in sorted(sequences.keys()):
	tick = time()
	with Pool(threads) as executor:
		executor.map(chromosome_crawl,[sequences[org][i] for i in sorted(sequences[org].keys())])
	print(f"Duration {time()-tick}")	
