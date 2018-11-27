#!/usr/bin/env python

import re
import os
import sys
import itertools

from collections import defaultdict
import argparse
from parsetrnas import *
from trnasequtils import *



#This program gets the mature tRNA sequences

parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--trnascan',  nargs='+', default=list(),
                   help='tRNAscan-SE file')
parser.add_argument('--rnacentral', nargs='+', default=list(),
                   help='RNAcentral tRNA file')
parser.add_argument('--genome', 
                   help='fasta sequence of genome')
parser.add_argument('--chromtranslate', 
                   help='translation table of chromosome names')
parser.add_argument('--bedfile', 
                   help='Output bedfile of coordinates for mature tRNA')
#parser.add_argument('--tag', nargs='1',
#                   help='tag to be added to tRNA name')
parser.add_argument('--maturetrnatable',
                   help='Output table of mature tRNAs')
parser.add_argument('--trnaalignment',
                   help='Output stockholm format alignment of mature tRNAs')


args = parser.parse_args()


#/projects/lowelab/users/holmes/pythonsource/trnastuff/ucscToINSDC.txt

'''
/projects/lowelab/users/holmes/pythonsource/trnaseq/getmaturetrnas.py --rnacentral hg19-trnacentral.txt  --genome /scratch/encodeshortrna/hg19.fa --chromtranslate /projects/lowelab/users/holmes/pythonsource/trnatest/trnastuff/ucscToINSDC.txt --maturetrnatable=hg19-trnatable.txt --trnaalignment hg19-trnaalign.stk --maturetrnatable=sacCer3-trnatable.txt

/projects/lowelab/users/holmes/pythonsource/trnaseq/getmaturetrnas.py --rnacentral /projects/lowelab/users/holmes/pythonsource/trnatest/hg19-trnacentral.txt  --genome /scratch/encodeshortrna/hg19.fa --chromtranslate /projects/lowelab/users/holmes/pythonsource/trnatest/trnastuff/ucscToINSDC.txt --maturetrnatable=hg19-trnatable.txt --trnaalignment=hg19-trnaalign.stk --maturetrnatable=sacCer3-trnatable.txt

Need to add bits to create maf files and bed files for these genomes

/projects/lowelab/users/holmes/pythonsource/trnaseq/getmaturetrnas.py --rnacentral sacCer3-trnacentral.txt  --genome sacCer3.fa --chromtranslate NameConversion.txt --bedfile=sacCer3-maturetRNAs.bed --maturetrnatable=sacCer3-trnatable.txt >sacCer3-maturetrnas.fa
'''



alltrnas = list()
trnascantrnas = list()
for currfile in args.trnascan:
    trnascantrnas.extend(readtRNAscan(currfile, args.genome))
trnacentraltrnas = list()
for currfile in args.rnacentral:
    trnacentraltrnas.extend(readrnacentral(currfile,args.chromtranslate,mode = 'transcript'))
    
alltrnas = list(getuniquetRNAs(trnascantrnas)) + trnacentraltrnas

trnabed = None
if args.bedfile:
    trnabed = open(args.bedfile, "w")

trnatable = None
if args.maturetrnatable:
    trnatable = open(args.maturetrnatable, "w")




def readmultistk(struct):
    currrecord = ""
    structs = list()
    for line in struct.split("\n"):
        currrecord += line+"\n"
        if line == "//":
            yield currrecord
            currrecord = ""
            

margin = 20
anticodoncount = defaultdict(int)
trnanames = dict()
trnalist = list()
for currtrans in alltrnas:
    if currtrans.name is None:
        name = 'tRNA-'+currtrans.amino + currtrans.anticodon+ str(anticodoncount[currtrans.anticodon]+ 1)
        currtrans.name = name
        anticodoncount[currtrans.anticodon] += 1
        
        
scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"            
trnacmfile = scriptdir+'trnamature-euk.cm'
stkfile = args.trnaalignment
if args.trnaalignment:
    devnull = open(os.devnull, 'w')
    seqfile = tempmultifasta(((currtrans.name, currtrans.getmatureseq()) for currtrans in alltrnas))
    cmcommand = ['cmalign', "-o", stkfile,"--nonbanded", "-g",trnacmfile,seqfile.name]
    #print >>sys.stderr, " ".join(cmcommand)
    cmrun = subprocess.Popen(cmcommand, stdout = devnull)
    result = cmrun.wait()
    if result:
        print >>sys.stderr, "Failure to align tRNAs"
        sys.exit(1)
    #stkout = cmrun.communicate()[0]
    #trnaalign = readrnastk(stkout.split("\n"))[0]
    seqfile.close()
    
    
for currtrans in alltrnas:
    name = currtrans.name
    #trnanames[name] = currtrans
    trnalist.append(name)
    #print >>sys.stderr, name
    print ">"+name
    print str("N" * margin) +currtrans.getmatureseq()+str("N" * margin)
    if trnatable is not None:
        print >>trnatable, "\t".join([name,",".join(currlocus.name for currlocus in currtrans.loci),currtrans.amino,currtrans.anticodon])
    if trnabed is not None:
        transcriptrange = GenomeRange("genome", name, margin, margin + len(currtrans.getmatureseq()), strand = "+", name = name)
        print >>trnabed, transcriptrange.bedstring()
    


#sys.exit()
trnamods = dict()
allmods = set()
nomatches = 0
trnamismatches = dict()

