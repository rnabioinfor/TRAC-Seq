#!/usr/bin/env python

import pysam
import os
import sys
import argparse
import subprocess
from collections import defaultdict
from trnasequtils import *


parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--trnaloci',
                   help='bed file of tRNA loci')      
parser.add_argument('--genomefile',
                   help='fasta file of genome')
parser.add_argument('--stkfile',
                   help='stockholm output file')


'''
~/pythonsource/aligntrnalocus.py --genomefile= --stkfile= --trnaloci

~/pythonsource/trnaseq/aligntrnalocus.py --genomefile=hg19.fa --stkfile=hg19-trnaloci.stk --trnaloci=hg19-trnaloci.bed
'''
args = parser.parse_args()

stkfile = args.stkfile
genomefile = args.genomefile


scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"            
trnacmfile = scriptdir+'TRNAinf-euk.cm'

trnaloci = list(readbed(args.trnaloci, orgdb = "genome", seqfile=genomefile))
lociseqs = getseqdict(trnaloci, faifiles = {"genome":genomefile+".fai"})
#print lociseqs
#lociseqfile = tempmultifasta(lociseqs)
devnull = open(os.devnull, 'w')
seqfile = tempmultifasta(lociseqs.iteritems())
cmcommand = ['cmalign', "-o", stkfile,"--nonbanded", "-g",trnacmfile,seqfile.name]
#print >>sys.stderr, " ".join(cmcommand)
cmrun = subprocess.Popen(cmcommand, stdout = devnull)
result = cmrun.wait()
if result:
    print >>sys.stderr, "Failure to align tRNAs"
    sys.exit(1)
seqfile.close()
devnull.close()

#trnaalign
    