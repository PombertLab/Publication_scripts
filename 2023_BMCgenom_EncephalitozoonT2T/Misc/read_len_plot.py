#!/usr/bin/python
## Pombert lab, 2022
version = '0.5d'
updated = '2022-07-05'
name = 'read_len_plot.py'

import os
import sys
import gzip
import argparse
import matplotlib.pyplot as plt
import re

################################################################################
## README
################################################################################

usage = f"""
NAME		{name}
VERSION		{version}
UPDATED		{updated}
SYNOPSIS	Plots the read length distribution for a given FASTQ dataset with 
		matplotlib

COMMAND		{name} \\
		  -f reads.fastq \\
		  -c darkorange \\
		  -o read_distribution.svg read_distribution.pdf \\
		  -x 50000

I/O OPTIONS:
-f (--fastq)	FASTQ file to plot (GZIP files are supported)
-d (--outdir)	Output directory [Default: ./]
-m (--metrics)	Metrics output file [Default: read_metrics.txt]
-o (--output)	Save plot to specified output file(s)
		## Defaults to matplotlib GUI otherwize
		## Supported formats: png, pdf, png, ps and svg

PLOT OPTIONS:
-c (--color)	Color to use; red, green, blue... [Default: green]
		# https://matplotlib.org/stable/gallery/color/named_colors.html
-b (--bar)	Bar type: Read sum or read count [Default: sum]
-h (--height)	Figure height in inches [Default: 10.8]
-w (--width)	Figure width in inches [Default: 19.2]
-x (--xmax)	Set max X-axis value [Default: automatic]
-t (--ticks)	Set ticks every X kb [Default: 5]
-y (--yscale)	Set yscale: linear or log [Default: linear]
--title		Set title; defaults to file basename if not set
--title_font	Set title font: normal, bold, heavy [Default: normal]
"""

# Print custom message if argv is empty
if (len(sys.argv) <= 1):
	print(usage)
	sys.exit()

################################################################################
## Create command lines switches
################################################################################

cmd = argparse.ArgumentParser(add_help=False)
cmd.add_argument("-f", "--fastq")
cmd.add_argument("-o", "--output", nargs='*')
cmd.add_argument("-d", "--outdir", default='./')
cmd.add_argument("-m", "--metrics", default='read_metrics.txt')
cmd.add_argument("-c", "--color", default='green')
cmd.add_argument("-b", "--bar", default='sum', choices=['sum', 'count'])
cmd.add_argument("-h", "--height", default=10.8)
cmd.add_argument("-w", "--width", default=19.2)
cmd.add_argument("-x", "--xmax", type=int)
cmd.add_argument("-t", "--ticks", type=int, default=5)
cmd.add_argument("-y", "--yscale", default='linear', choices=['linear', 'log'])
cmd.add_argument("--title")
cmd.add_argument("--title_font", default='normal', choices=['normal', 'bold', 'heavy'])
args = cmd.parse_args()

fastq = args.fastq
output = args.output
outdir = args.outdir
metrics_file = args.metrics
bar = args.bar
height = args.height
width = args.width
rgb = args.color
xmax = args.xmax
set_ticks = args.ticks
yscale = args.yscale
title = args.title
title_font = args.title_font


################################################################################
## Working on output directory
################################################################################

if output is not None:
	if os.path.isdir(outdir) == False:
		try:
			os.makedirs(outdir)
		except:
			sys.exit(f"Can't create directory {outdir}...")


################################################################################
## Working on FASTQ file
################################################################################

read_sizes = []
line_counter = 0

FH = None

## Check if file is gzipped or not
try:
	FH = gzip.open(fastq,'r')
except:
	FH = open(fastq,'r')
else:
	FH = gzip.open(fastq,'r')

print(f"\nWorking on {fastq}...\n")

for line in FH:
	line_counter += 1
	if (line_counter == 2):
		line = line.strip()
		read_size = len(line)
		read_sizes.append(read_size)
	elif (line_counter == 4):
		line_counter = 0

################################################################################
## Calculate read metrics
################################################################################

read_num = "{:,}".format(len(read_sizes))
read_sum = "{:,}".format(sum(read_sizes))
longest_read = "{:,}".format(max(read_sizes))
shortest_read = "{:,}".format(min(read_sizes))

average = sum(read_sizes)/len(read_sizes)
average = int(round(average))
average = "{:,}".format(average)

median_location = int(round(len(read_sizes)/2))
median = read_sizes[median_location]
median = "{:,}".format(median)

# Function to calculate n metrics; e.g n50, n75, n90
def n_metric(list, n):
	n_threshold = int(sum(list)*n)
	n_sum = 0
	nmetric = 0

	for x in list:
		n_sum += x
		if n_sum >= n_threshold:
			nmetric = x
			break

	nmetric = "{:,}".format(nmetric)
	return nmetric

# n50, 75, 90
read_sizes.sort(reverse=True)
n50 = n_metric(read_sizes,0.5)
n75 = n_metric(read_sizes,0.75)
n90 = n_metric(read_sizes,0.9)

# Print metrics to file or STDOUT
pmetrics = f"""Metrics for {fastq}:

Total bases:\t{read_sum}
# reads:\t{read_num}
Longest:\t{longest_read}
Shortest:\t{shortest_read}
Average:\t{average}
Median:\t\t{median}
N50:\t\t{n50}
N75:\t\t{n75}
N90:\t\t{n90}
"""

if output is None:
	print(pmetrics)
else:
	metrics_output = outdir + '/' + metrics_file
	METRICS = open(metrics_output,'w')
	print(pmetrics, file=METRICS)


################################################################################
## Plot read length distribution with matplotlib (using bar chart)
################################################################################

##### Metrics text box

adjust_l = 0
metrics_list = [read_sum, read_num, longest_read, shortest_read, average, median, n50, n75, n90]
for metric in metrics_list:
	if len(metric) > adjust_l:
		adjust_l = len(metric)

metrics = f"""
	Total bases: {read_sum.rjust(adjust_l + 1)}
	# reads: {read_num.rjust(adjust_l + 1)}
	Longest: {longest_read.rjust(adjust_l + 1)}
	Shortest: {shortest_read.rjust(adjust_l + 1)}
	Average: {average.rjust(adjust_l + 1)}
	Median: {median.rjust(adjust_l + 1)}
	N50: {n50.rjust(adjust_l + 1)}
	N75: {n75.rjust(adjust_l + 1)}
	N90: {n90.rjust(adjust_l + 1)}
"""
metrics = metrics.replace("\t","")

##### Bins, ticks and text box location

# Bin size + dictionary to store read size distributions
binsize = 1000
reads_distr = {
	'sum': {},
	'count': {}
}

# Checking for max value, either from data or the command line
max_val = None
if xmax is None:
	max_val = max(read_sizes)
else:
	max_val = xmax

# Creating bins
num_bins = int(max_val/binsize) + 1
key_labels = []

# Creating bin labels
for x in range(0,num_bins):
	label = f"{x}k"
	key_labels.append(label)
	reads_distr['sum'][label] = 0
	reads_distr['count'][label] = 0

# Calculating sum of all bases (in Mb) + read count per bin
for size in read_sizes:
	bin_loc = int(size/binsize)
	bin_loc = f"{bin_loc}k"
	read_mb = size/1000000
	if bin_loc in reads_distr[bar].keys():
		reads_distr['sum'][bin_loc] += read_mb
		reads_distr['count'][bin_loc] += 1
	else:
		reads_distr['sum'][bin_loc] = read_mb
		reads_distr['count'][bin_loc] = 1

# Setting ticks for plot
ticks = []
labels = []
for i in range(0, num_bins, set_ticks):
	ticks.append(i)
	labels.append(key_labels[i])

# Metrics text box location (top right corner)
x_metrics_location = ticks[-1] - 1
max_bin_value = 0
for key in reads_distr[bar].keys():
	if reads_distr[bar][key] > max_bin_value:
		max_bin_value = reads_distr[bar][key]
y_metrics_location = max_bin_value - 1

##### Plotting bar chart

# Setting default image to widescreen by default
plt.rcParams["figure.figsize"] = (width,height)

# Setting bar labels
xlabel = 'Read sizes'
ylabel = None
if bar == 'sum':
	ylabel = 'Total bases (in Mb)'
elif bar == 'count':
	ylabel = 'Total read count'

if title is None:
	title = os.path.basename(fastq)

plt.title(
	title,
	loc='center',
	fontsize = 12,
	y = 1.0,
	pad=-50,
	fontweight=title_font
)
plt.text(
	x_metrics_location,
	y_metrics_location,
	metrics,
	fontsize=10,
	family='monospace',
	va='top',
	ha='right'
)

plt.xlabel(xlabel)
plt.xlim(0,num_bins)
plt.xticks(ticks,labels)
plt.ylabel(ylabel)
plt.yscale(yscale)
plt.bar(list(reads_distr[bar].keys()), reads_distr[bar].values(), color=rgb, align='edge')

# Output either to matplotlib GUI or file
if output is None:
	plt.show()
else:
	for x in output:
		filename = outdir + '/' + x
		print(f"  Creating {filename}...")
		plt.savefig(filename)
