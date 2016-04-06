#!/usr/bin/env python
import sys
import subprocess
import pysam

if len(sys.argv) < 3:
    sys.stderr.write('usage:\t' + \
                     sys.argv[0] + \
                     ' <ped> <padding (e.g. +-500)> <sample ids (opt)>\n')
    sys.exit(1)

ped_file_name = sys.argv[1]
padding = int(sys.argv[2])
samples = []
if len(sys.argv) == 4:
    samples = sys.argv[3:]

f = open(ped_file_name, 'r')

sample_to_bam = {}
for l in f:
    A = l.rstrip().split('\t')
    if A[0] == 'family_id':
        continue
    if (len(samples) == 0)  or (A[1] in samples):
        sample_to_bam[A[1]] = [A[-1],0]
f.close()

infile = pysam.VariantFile("-", "r")

i = 0
for sample in infile.header.samples:
    if infile.header.samples[i] in sample_to_bam.keys():
        sample_to_bam[infile.header.samples[i]][1] = i
    i += 1

for variant in infile:
    if (variant.info['SVTYPE'] != 'DEL')  and (variant.info['SVTYPE'] != 'DUP'):
        continue

    print 'var_' + variant.id + '.txt'
    f = open('var_' + variant.id + '.txt', 'w')
    #print variant.info['SVTYPE'], variant.chrom, variant.start, variant.stop
    f.write(' '.join([str(x) for x in[variant.info['SVTYPE'],\
                                      variant.chrom, \
                                      max(0,variant.start-padding), \
                                      variant.stop + padding]]) + '\n')
    f.write(str(padding) + '\n')
    for sample in sample_to_bam:
        gt_raw = variant.samples[sample_to_bam[sample][1]].allele_indices
        gt = '/'.join([str(x) for x in gt_raw])
        cmd = "samtools depth -r " + \
               variant.chrom + ":" + \
                    str(max(0,variant.start-padding)) + "-" + \
                    str(variant.stop+padding) + " " + \
              sample_to_bam[sample][0] + \
              "| cut -f2,3"
              #"| cut -f2,4"
              #"| tr '\\n' ' '"
        #cmd = ["ls", "-l", sample_to_bam[sample][0] ]
        #print cmd
        p = subprocess.Popen( [cmd], stdout=subprocess.PIPE,shell=True )
        depth=[sample, str(gt)]
        last  = 0
        for line in p.stdout:
            pos,dep = [int(x) for x in line.rstrip().split()]
            if (last != 0) and (pos != last +1):
                for z in range(0,pos-last-1):
                    depth.append('0')
            else:
                depth.append(str(dep))
            last = pos
        f.write(' '.join(depth) + '\n')
    f.close()
