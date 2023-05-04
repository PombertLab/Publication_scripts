#!/usr/bin/python3

name = 'k_plotter.py'
version = '0.2.1'
updated = '2022-11-22'

usage = f'{name} -i kmer_redundancy -o outdir -p gene_positions'

import argparse
from re import match
from sys import argv
from os import listdir,makedirs
from os.path import isdir,exists
from matplotlib import pyplot as plt

if len(argv) < 2:
	print(f"\n{usage}\n\n")
	exit()

GetOptions = argparse.ArgumentParser()

GetOptions.add_argument('-i','--indir')
GetOptions.add_argument('-o','--outdir',default = "KmerFigures")
GetOptions.add_argument('-p','--positions',nargs="*")

GetOptions = GetOptions.parse_args()

indir = GetOptions.indir
outdir = GetOptions.outdir
positions = GetOptions.positions

if not exists(outdir):
	makedirs(outdir)

for org in listdir(indir):

	if isdir(f"{indir}/{org}"):
		print(f"Organism: {org}")
		for chromosome in sorted(listdir(f"{indir}/{org}")):
			if isdir(f"{indir}/{org}/{chromosome}"):
				if not exists(f"{outdir}/{chromosome}"):
					makedirs(f"{outdir}/{chromosome}")
				print(f"Chromosome: {chromosome}")
				data = {'gc':{'x':[],'y':[]},'kmer':{},'type':{'rRNA':[],'KNOWN':[],'HYPO':[]}}
				for file in sorted(listdir(f"{indir}/{org}/{chromosome}")):
					if ".kmer_count" in file:
						kmer_size, window_size, slide = [int(i) for i in match(".*\.kmer_(\d+).window_(\d+).slide_(\d+).kmer_count",file).groups()]
						print(f"Kmer Size: {kmer_size}")
						IN = open(f"{indir}/{org}/{chromosome}/{file}","r")
						data['kmer'][kmer_size] = {'window_size':window_size,'slide':slide,'x':[],'rbuk':[],'rbk50':[],'rbaak':[],'ar':[]}
						for line in IN:
							line = line.replace("\n","")
							if line[0] != '#':
								values = line.split("\t")
								x,rbuk,rbk50,rbaak,ar = [float(i) for i in values]
								data['kmer'][kmer_size]['x'].append(x)
								data['kmer'][kmer_size]['rbuk'].append(rbuk)
								data['kmer'][kmer_size]['rbk50'].append(rbk50)
								data['kmer'][kmer_size]['rbaak'].append(rbaak)
								data['kmer'][kmer_size]['ar'].append(ar)
						IN.close()

					if ".gc_content" in file:
						IN = open(f"{indir}/{org}/{chromosome}/{file}","r")
						for line in IN:
							line = line.replace("\n","")
							if line[0] != '#':
								values = line.split("\t")
								x,c = [int(values[0]),float(values[1])]
								data['gc']['x'].append(x)
								data['gc']['y'].append(c)
					
				if f"{chromosome}.gene_loc" in positions:
					IN = open(f"{chromosome}.gene_loc","r")
					for line in IN:
						line = line.replace("\n","")
						type,start,end = line.split("\t")
						data['type'][type].append([int(start),int(end)])

				for kmer_size in sorted(data['kmer'].keys()):
					plt.figure(figsize=(16,9))
					ax = plt.subplot()
					d1 = ax.plot(data['gc']['x'],[i*100 for i in data['gc']['y']],label="GC Content",color='orange') #color="#03658c",
					ax.set_ylabel("% GC content")
					ax.set_ylim([0,100])
					
					for type in data['type']:
						for start,end in data['type'][type]:
							color = ''
							if type == 'rRNA':
								color = 'purple'
							if type == 'KNOWN':
								color = 'lightseagreen'
							if type == 'HYPO':
								color = 'steelblue'
							plt.fill_between([start,end],0.0,2.5,color=color,alpha=0.25)

					ax2 = ax.twinx()
					d2 = ax2.plot(data['kmer'][kmer_size]['x'],data['kmer'][kmer_size]['rbuk'],label=f"kmer (size={kmer_size}) Redundancy") #color='#730220',
					# ax2.plot(data['gc']['x'],[0 for i in data['gc']['y']],color="grey",linestyle="--",linewidth=1)
					ax2.set_ylim([0,1])
					ax2.set_xlim([data['kmer'][kmer_size]['x'][0],data['kmer'][kmer_size]['x'][-1]])
					ax2.set_xlabel("Chromosome Position (bps)")
					ax2.set_ylabel("Proprotion of Redundancy")
					Ins = d1 + d2
					labels = [l.get_label() for l in Ins]
					plt.legend(Ins,labels)
					plt.title(f"{org}: {chromosome}")
					plt.savefig(f"{outdir}/{chromosome}/{chromosome}.kmer_{kmer_size}.window_{window_size}.slide_{slide}.png",dpi=300,format='png')
					plt.savefig(f"{outdir}/{chromosome}/{chromosome}.kmer_{kmer_size}.window_{window_size}.slide_{slide}.svg",dpi=300,format='svg')
					plt.close()