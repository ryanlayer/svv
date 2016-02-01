#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
import pylab
import random
import matplotlib.gridspec as gridspec
from optparse import OptionParser
from operator import itemgetter
import subprocess

delim = '\t'
parser = OptionParser()

parser.add_option("-l",
                  "--log_y",
                  action="store_true", dest="logy", default=False,
                  help="Use log scale for y-axis")

parser.add_option("-o",
                  "--output_file",
                  dest="output_file",
                  help="Data file")

parser.add_option("--y_min",
                  dest="min_y",
                  help="Min y value")

parser.add_option("--y_max",
                  dest="max_y",
                  help="Max y value")

parser.add_option("--line_style",
                  dest="line_style",
                  default=",",
                  help="Line style")

parser.add_option("-e",
                  "--exon_file",
                  dest="exon_file",
                  help="Exon file")

parser.add_option("-n",
                  "--names",
                  dest="names",
                  help="Gene names")



(options, args) = parser.parse_args()
if not options.output_file:
    parser.error('Output file not given')

trans_gene_names = {}
if options.names:
    f = open(options.names, 'r')
    for l in f:
        A=l.rstrip().split()
        trans_gene_names[A[0]] = A[1]

names = []
gt = []
Y=[]
title=''
padding=''
for l in sys.stdin:
    if title == '':
        title = l.rstrip()
        continue

    if padding == '':
        padding = int(l.rstrip())
        continue

    A = l.split()
    names.append(A[0])
    gt.append(A[1])
    Y.append([float(x) for x in A[2:]])


(chrm,x_start,x_end) = title.split()[1:]

cmd = "tabix " + options.exon_file + " " + \
        chrm +":"+x_start+"-"+x_end
print cmd
p = subprocess.Popen( [cmd], stdout=subprocess.PIPE,shell=True )
transcripts = {}
for line in p.stdout:
    A = line.rstrip().split('\t')
    (transcript_name, foo, exon_i, bla, chrom, start, strand) = A[3].split('_')
    if transcript_name not in transcripts:
        transcripts[transcript_name] = []

    print int(A[1]), int(x_start)
    transcripts[transcript_name].append([ int(exon_i),
                                          max(0,int(A[1]) - int(x_start)),
                                          max(0,int(A[2]) - int(x_start))])

print transcripts

matplotlib.rcParams.update({'font.size': 12})
fig = matplotlib.pyplot.figure(figsize=(5,10),dpi=300)
fig.subplots_adjust(wspace=.05,left=.01,bottom=.01)

gs = gridspec.GridSpec(len(Y) + 1, 1, height_ratios = [1] + [2] * len(Y))

ax = matplotlib.pyplot.subplot(gs[0])

i = 1
for transcript in transcripts:
    min_p = transcripts[transcript][0][1]
    max_p = transcripts[transcript][-1][2]
    ax.plot([min_p,max_p],[i/10.0,i/10.0], color='black')
    for exon in sorted(transcripts[transcript], key=itemgetter(0)):
        ax.plot([exon[1],exon[2]],[i/10.0,i/10.0], color='black', lw=3)

    out_name = transcript
    if options.names and transcript in trans_gene_names:
        out_name = trans_gene_names[transcript]
    ax.text(min_p,i/10.0,out_name,fontsize=8)
    i+=1
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.xaxis.set_ticks_position('none') 
ax.yaxis.set_ticks_position('none') 


y_max = max(max(y) for y in Y)
x_len = len(Y[0])
N = len(Y)

i = 1

fig.suptitle(title)
for y in Y:
    #ax = fig.add_subplot(N,1,i)
        #ax.plot([100, 300], [max(y), max(y)], 'k-', lw=5)
    ax = matplotlib.pyplot.subplot(gs[i])
    ax.text(0,y_max,names[i-1] + " " + gt[i-1])
    ax.plot([padding,padding],[y_max,0], color='black')
    ax.plot([len(y) - padding,len(y) - padding],[y_max,0], color='black')

    ax.plot(range(len(y)),y,options.line_style,color='black', linewidth=1)

    if options.logy:
        ax.set_yscale('log')

    #ax.set_yticklabels([])
    ax.set_xticklabels([])

    if ((options.max_y) and (options.min_y)):
        ax.set_ylim(float(options.min_y),float(options.max_y))

    ax.set_ylim(float(0),float(y_max))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


    i+=1

matplotlib.pyplot.savefig(options.output_file,bbox_inches='tight')
