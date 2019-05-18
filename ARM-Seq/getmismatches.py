#!/usr/bin/env python

import pysam
import sys
import argparse
import string
import itertools
from collections import defaultdict
from trnasequtils import *

dbname = 'sacCer3'
'''
~/pythonsource/trnaseq/getmismatches.py --samplefile=YeastAging.txt --bedfile=sacCer3-maturetRNAs.bed --stkfile=sacCer3-trnaalign.stk >aging-mismatches.txt
'''

gapchars = set("-._~")

count = 0

class readmismatch:
    def __init__(self, region):
        self.mismatches = list()
                
    def coveragelist(self):
        return self.mismatches
    def coveragealign(self, alignment, gapoutput = "NA",sizefactor = 1):
        if len(self.mismatches) != len(string.translate(alignment, None, str(gapchars))):
            print >>sys.stderr, "Alignment length does not match bed length"            
        i = 0
        for curr in alignment:
            #print >>sys.stderr, curr
            if curr in gapchars:
                yield gapoutput
            else:
                yield self.mismatches[i]/sizefactor
                i += 1



parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--bedfile',  nargs='+', default=list(),
                   help='bed file with mature tRNA features')
parser.add_argument('--samplefile',
                   help='Sample file in format')
parser.add_argument('--sizefactors',
                   help='Optional file including size factors that will be used for normalization')
parser.add_argument('--stkfile',
                   help='Stockholm file')




args = parser.parse_args()

sampledata = samplefile(args.samplefile)
samples = sampledata.getsamples()

alltrnas = list()
#print >>sys.stderr, " ".join(bamlist)
trnastk = list(readrnastk(open(args.stkfile, "r")))[0]
try:
    trnalist = list()                       
    for currfile in args.bedfile:
        trnalist.extend(list(readbed(currfile)))       
except IOError as e:
    print >>sys.stderr, e
    sys.exit()


#featurelist = list(readbed(sys.argv[1]))
#featurelist = list(readbed(sys.argv[1]))

#./countcomplete.py hg19-nontrnas.bed hg19-tRNAs.bed hg19-complete-tRNApad.fa

sizefactor = defaultdict(lambda: 1)
if args.sizefactors:
    sizefactor = getsizefactors(args.sizefactors)
    
    
'''
./getcoverage.py ../combinedb/sacCer3-fatRNAs.bed sacCer3-agingtranscripts.bed >sacCer3-agingcount.txt
'''

featcount = defaultdict(int)
#featurelist = list(curr for curr in featurelist if curr.name == 'unknown20')
#print >>sys.stderr, "***"
#lengths = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))


maxoffset = 10
allcoverages = dict()


def nasum(operands, naval = "NA"):
    if sum("NA"== curr for curr in operands) == len(operands):
        return "NA"
    elif not any("NA"== curr for curr in operands):
        return sum(operands)
    else:
        print >>sys.stderr, "Trying to add incompatible alignments"
        sys.exit(1)
    #return ",".join(str(curr) for curr in operands)
def sumsamples(coverage,sampledata, repname, currfeat, sizefactors = defaultdict(lambda: 1)):
    return (nasum(curr) for curr in itertools.izip(*(allcoverages[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name],sizefactor = sizefactor[currsample]) for currsample in sampledata.getrepsamples(repname))))

    
    
for currsample in samples:
    currbam = sampledata.getbam(currsample)
    allcoverages[currsample] = dict()
    try:
        #print >>sys.stderr, currbam
        pysam.index(""+currbam)
        bamfile = pysam.Samfile(""+currbam, "rb" )  
    except IOError as ( strerror):
        print >>sys.stderr, strerror
        sys.exit()
        
    for currfeat in trnalist:
        allcoverages[currsample][currfeat.name] = readmismatch(currfeat)
        for reference, basecounts in getpileuprange(bamfile, currfeat):
            matchcount = basecounts[reference]
            
            mismatchcount = sum(basecounts.values()) - matchcount
            allcoverages[currsample][currfeat.name].mismatches.append(mismatchcount)
            
for currfeat in trnalist:
    for currsample in samples:
        print currfeat.name+"\t"+currsample+"\t"+"\t".join(str(curr) for curr in allcoverages[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name],sizefactor = sizefactor[currsample]))
        