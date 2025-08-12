#!/usr/bin/env python3
# Olivier CREPEAULT

name = 'plot_Coverage.py'
version = '0.2.0'
updated = '2025-08-12'

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

################################################################################
## README
################################################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Plots RNASeq coverage for a specified chromosome.

REQS        matplotlib

OPTIONS:
-f (--fasta)        Chromosome sequence in FASTA format
-c (--coverage)     RNASeq coverage file from samtools depth -aa
-w (--window)       Sliding window size [Default: 10000]
-s (--step)         Step size for sliding window [Default: 1000]
-o (--output)       Output directory to save the plot [Default: ./]
-y (--ymax)         Set a maximum y-axis value for uniformity
--format            Output file format (e.g., png, svg) [Default: png]
--version           Show script version
"""

# Print custom message if argv is empty
if (len(sys.argv) <= 1):
    print(usage)
    exit(0)

################################################################################
## Create command lines switches
################################################################################
    
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-f', '--fasta')
parser.add_argument('-c', '--coverage')
parser.add_argument('-w', '--window', type=int, default=10000)
parser.add_argument('-s', '--step', type=int, default=1000)
parser.add_argument('-o', '--output', default='./')
parser.add_argument('-y', '--ymax', type=float)
parser.add_argument('--format', type=str, default='png')
parser.add_argument("--version", action='store_true')

args = parser.parse_args()

fasta = args.fasta
coverage = args.coverage
window = args.window
step = args.step
output = args.output
ymax = args.ymax
format = args.format
scversion = args.version

#########################################################################
### Version
#########################################################################

if scversion:
    print ("")
    print (f"Script:     {name}")
    print (f"Version:    {version}")
    print (f"Updated:    {updated}\n")
    exit(0)


################################################################################
## Functions
################################################################################

def read_rna_seq(rna_seq_coverage):
    
    """
    Function to read the RNASeq coverage file and extract the positions and coverage
    """

    rna_seq_df = pd.read_csv(rna_seq_coverage, sep='\t', header=None, names=['chromosome', 'position', 'coverage'])
    return rna_seq_df


def calculate_coverage(rna_seq_df, window_size, step_size):

    """
    Calculate coverage using sliding window
    """

    coverage = []
    total_length = rna_seq_df['position'].max()
    
    for i in range(1, total_length, step_size):
        window_data = rna_seq_df[(rna_seq_df['position'] >= i) & (rna_seq_df['position'] < i + window_size)]
        avg_coverage = window_data['coverage'].mean() if not window_data.empty else 0
        coverage.append(avg_coverage)
    
    return coverage


def save_plot(coverage, window_size, output_dir, fasta_file, step_size, y_max=None, file_format='png'):
    
    """
    Function to save the plot with standardized y-axis if provided
    """

    # Generate x-axis (window positions)
    x_axis_cov = np.arange(0, len(coverage)) * step_size

    # Plotting the data
    plt.figure(figsize=(10, 6))

    # Plot RNASeq coverage
    plt.plot(x_axis_cov, coverage, 'b-', label="Coverage")
    plt.xlabel('Position in chromosome')
    plt.ylabel('Coverage', color='b')
    
    # Set the y-axis limit if y_max is provided
    if y_max:
        plt.ylim(0, y_max)

    plt.title(f'RNASeq Coverage (Window: {window_size}, Step: {step_size})')

    # Get the base name of the fasta file and save the plot with sliding window and step size in the filename
    fasta_basename = os.path.splitext(os.path.basename(fasta_file))[0]
    output_path = os.path.join(output_dir, f"{fasta_basename}_coverage_window{window_size}_step{step_size}.{file_format}")
    
    # Save as specified format (png or svg, etc.)
    plt.savefig(output_path, format=file_format)
    print(f"Plot saved at: {output_path}")
    plt.close()


################################################################################
## Main Function
################################################################################

if __name__ == "__main__":

    # Create the output directory if it doesn't exist
    if not os.path.exists(output):
        os.makedirs(output)

    # Load RNASeq data
    rna_seq_data = read_rna_seq(coverage)

    # Calculate coverage using the step size
    coverage = calculate_coverage(rna_seq_data, window, step)

    # Save the plot with the sliding window, step size, and standardized y-axis if provided
    save_plot(coverage, window, output, fasta, step, ymax, format)
