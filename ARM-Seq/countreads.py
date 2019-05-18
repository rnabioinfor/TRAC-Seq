#!/usr/bin/env python

import pysam
import sys
import argparse
from collections import defaultdict
from trnasequtils import *

dbname = 'sacCer3'
'''
~/pythonsource/trnaseq/countreads.py --samplefile=/projects/lowelab/users/holmes/pythonsource/trnatest/KitComparison.txt --bedfile=hg19-transcripts.bed --maturetrnas=hg19-maturetRNAs.bed --trnaloci=hg19-trnas.bed

~/pythonsource/trnaseq/countreads.py --samplefile=/projects/lowelab/users/holmes/pythonsource/trnatest/testcomp.txt --bedfile=hg19-transcripts.bed --maturetrnas=hg19-maturetRNAs.bed --trnaloci=hg19-trnas.bed >allcounts.txt

Rscript ~/pythonsource/trnaseq/analyzecounts.R YeastAging featurecounts.txt agingshort.txt dmStat_Amino:dmAll_Amino dmStat_Amino:dmMet_Amino dmStat_Amino:dmLeu_Amino


This currently swaps when given huge bed files for features
'''

def getdupes(namelist):
    allset = set()
    for currname in namelist:
        if currname in allset:
            yield currname
        else:
            allset.add(currname)
        

count = 0



parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--samplefile',
                   help='Sample file in format')
parser.add_argument('--bedfile',  nargs='+', default=list(),
                   help='bed file with non-tRNA features')
parser.add_argument('--gtffile',  nargs='+', default=list(),
                   help='gtf file with non-tRNA features')
parser.add_argument('--gtftrnas',  nargs='+', default=list(),
                   help='gtf file with tRNA features')
parser.add_argument('--trnaloci',  nargs='+', default=list(),
                   help='bed file with tRNA features')
parser.add_argument('--maturetrnas',  nargs='+', default=list(),
                   help='bed file with mature tRNA features')



args = parser.parse_args()




wholetrnas = dict()
fivefrags = dict()
threefrags = dict()
trailerfrags = dict()
otherfrags = dict()
allfrags = dict()


alltrnas = list()


samplefiles = dict()

sampledata = samplefile(args.samplefile)
samples = sampledata.getsamples()



#print >>sys.stderr, " ".join(bamlist)
try:
    featurelist = list()
    trnaloci = list()
    #for currfile in args.bedfile:
    #    featurelist.extend(list(readfeatures(currfile)))
    trnalist = list()
    for currfile in args.trnaloci:
        trnaloci.extend(list(readbed(currfile)))
    for currfile in args.maturetrnas:
        trnalist.extend(list(readbed(currfile)))       
except IOError as e:
    print >>sys.stderr, e
    sys.exit()

#if checkduplicates(curr.name for curr in trnalist+trnaloci+featurelist):
    #sys.exit()
#featurelist = list(readbed(sys.argv[1]))
#featurelist = list(readbed(sys.argv[1]))

#./countcomplete.py hg19-nontrnas.bed hg19-tRNAs.bed hg19-complete-tRNApad.fa
'''
./countcomplete.py ../combinedb/sacCer3-fatRNAs.bed sacCer3-agingtranscripts.bed >sacCer3-agingcount.txt
'''

featcount = defaultdict(int)
allfeats = featurelist+trnaloci+trnalist

if len(set(curr.name for curr in allfeats)) < len(list(curr.name for curr in allfeats )):
    #print >>sys.stderr, list(curr.name for curr in featurelist )
    #print >>sys.stderr, len(set(curr.name for curr in featurelist))
    print >>sys.stderr, "Duplicate names in feature list:"
    #print >>sys.stderr, ",".join(getdupes(curr.name for curr in allfeats))
    #currname
    #sys.exit(1)


#featurelist = list(curr for curr in featurelist if curr.name == 'unknown20')
alltrnas = list(curr.name for curr in featurelist)
#print >>sys.stderr, "***"

#lengths = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
counts = defaultdict(lambda: defaultdict(int))
trnacounts = defaultdict(lambda: defaultdict(int))
trnawholecounts = defaultdict(lambda: defaultdict(int))
trnafivecounts = defaultdict(lambda: defaultdict(int))
trnathreecounts = defaultdict(lambda: defaultdict(int))

trnalocuscounts = defaultdict(lambda: defaultdict(int))
trnatrailercounts = defaultdict(lambda: defaultdict(int))

maxoffset = 10



for currsample in samples:
    currbam = sampledata.getbam(currsample)
    try:
        #print >>sys.stderr, currbam
        pysam.index(""+currbam)
        bamfile = pysam.Samfile(""+currbam, "rb" )  
    except IOError as ( strerror):
        print >>sys.stderr, strerror
        sys.exit()
    
    
    for currfeat in featurelist:
        try:
            for currread in getbamrange(bamfile, currfeat):
                if currfeat.coverage(currread) > 10:
                    counts[currsample][currfeat.name] += 1
        except ValueError:
            pass
    for currfeat in trnaloci:
        for currread in getbamrange(bamfile, currfeat):
            if currfeat.coverage(currread) > 10:
                trnalocuscounts[currsample][currfeat.name] += 1
    
    for currfeat in trnalist:
        if currfeat.name != "tRNA-Pro-CGG-1":
            pass
            #continue
        #print(currfeat)
        for currread in getbamrange(bamfile, currfeat):
            
            if not currfeat.strand == currread.strand:
                continue
            if not currfeat.coverage(currread) > 10:
                continue
                #print >>sys.stderr, currsample
                #print >>sys.stderr, currread.bedstring()
                #print >>sys.stderr, currfeat.bedstring()
                #print >>sys.stderr, "********"
                
            trnacounts[currsample][currfeat.name] += 1
                
            fragtype = getfragtype(currfeat, currread)
            if fragtype == "Whole":
                trnawholecounts[currsample][currfeat.name] += 1
            elif fragtype == "Fiveprime":
                trnafivecounts[currsample][currfeat.name] += 1
            elif fragtype == "Threeprime":
                trnathreecounts[currsample][currfeat.name] += 1
              
                    #print >>sys.stderr, str(currread.start - currfeat.start)+"-"+str(currread.end - currfeat.start)  
                    #print >>sys.stderr, str(currfeat.start - currfeat.start)+"-"+str(currfeat.end - currfeat.start)
                    #print >>sys.stderr, "****"
                        



print "\t".join(samples)

trnanames = set()
for currfeat in trnalist:

    print currfeat.name+"\t"+"\t".join(str(trnacounts[currsample][currfeat.name]) for currsample in samples)
    #print currfeat.name+"_wholecounts\t"+"\t".join(str(trnawholecounts[currsample][currfeat.name]) for currsample in samples)
    #print currfeat.name+"_fiveprime\t"+"\t".join(str(trnafivecounts[currsample][currfeat.name]) for currsample in samples)
    #print currfeat.name+"_threeprime\t"+"\t".join(str(trnathreecounts[currsample][currfeat.name]) for currsample in samples)
    #print currfeat.name+"_other\t"+"\t".join(str(trnacounts[currsample][currfeat.name] - (trnathreecounts[currsample][currfeat.name] + trnafivecounts[currsample][currfeat.name] + trnawholecounts[currsample][currfeat.name])) for currsample in samples)
    
    

for currfeat in trnaloci:
    print currfeat.name+"\t"+"\t".join(str(trnalocuscounts[currsample][currfeat.name]) for currsample in samples)
    
for currfeat in featurelist:
    if currfeat.name in trnanames:
        continue
    trnanames.add(currfeat.name)
    print currfeat.name+"\t"+"\t".join(str(counts[currsample][currfeat.name]) for currsample in samples)


    
    
