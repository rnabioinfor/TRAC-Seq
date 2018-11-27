#!/usr/bin/env python

import pysam
import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *

dbname = 'sacCer3'
'''


~/pythonsource/trnaseq/countreadtypes.py --samplefile=agingshort.txt   --maturetrnas=sacCer3-maturetRNAs.bed --trnaloci=sacCer3-trnas.bed --bedfile sacCer3-rRNA.bed sacCer3-snoRNAs.bed sacCer3-transcripts.bed 


~/pythonsource/trnaseq/countreadtypes.py --sizefactors=tRNAseq_SizeFactors.txt --combinereps --samplefile=agingshort.txt   --maturetrnas=sacCer3-maturetRNAs.bed --countfrags --trnaloci=sacCer3-trnas.bed --bedfile sacCer3-rRNA.bed sacCer3-snoRNAs.bed  >repfragtypenormcounts.txt
'''



count = 0



parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--samplefile',
                   help='Sample file in format')
parser.add_argument('--sizefactors',
                   help='Optional file including size factors that will be used for normalization')
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
parser.add_argument('--countfrags', action="store_true", default=False,
                   help='Seperate tRNA fragment types')
parser.add_argument('--combinereps', action="store_true", default=False,
                   help='Sum samples that are replicates')

parser.add_argument('--trnanormfile',
                   help='Create normalization file to use to normalize to total tRNA reads')
parser.add_argument('--allreadsnormfile',
                   help='Create normalization file to use to normalize to total reads')


args = parser.parse_args()


countfrags = args.countfrags
combinereps = args.combinereps

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

sizefactor = defaultdict(lambda: 1)
if args.sizefactors:
    sizefactor = getsizefactors(args.sizefactors)
    

#print >>sys.stderr, " ".join(bamlist)
try:
    featurelist = dict()
    trnaloci = dict()
    trnalist = dict()
    for currfile in args.bedfile:
        featurelist[currfile] = RangeBin(readfeatures(currfile))
    
    for currfile in args.trnaloci:
        trnaloci[currfile] = RangeBin(readbed(currfile))
    for currfile in args.maturetrnas:
        trnalist[currfile] = RangeBin(readbed(currfile))
except IOError as e:
    print >>sys.stderr, e
    sys.exit()


#featurelist = list(readbed(sys.argv[1]))
#featurelist = list(readbed(sys.argv[1]))

#./countcomplete.py hg19-nontrnas.bed hg19-tRNAs.bed hg19-complete-tRNApad.fa
'''
./countcomplete.py ../combinedb/sacCer3-fatRNAs.bed sacCer3-agingtranscripts.bed >sacCer3-agingcount.txt
'''

featcount = defaultdict(int)




#featurelist = list(curr for curr in featurelist if curr.name == 'unknown20')
#print >>sys.stderr, "***"

#lengths = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
counts = defaultdict(lambda: defaultdict(int))
trnacounts = defaultdict(lambda: defaultdict(int))
trnawholecounts = defaultdict(lambda: defaultdict(int))
trnafivecounts = defaultdict(lambda: defaultdict(int))
trnathreecounts = defaultdict(lambda: defaultdict(int))
trnatrailercounts  = defaultdict(lambda: defaultdict(int))

trnalocuscounts = defaultdict(lambda: defaultdict(int))
trnatrailercounts = defaultdict(lambda: defaultdict(int))
othercounts = defaultdict(int)

maxoffset = 10

bedlist = list(featurelist.iterkeys())
locilist = list(trnaloci.iterkeys())
'''
I use readbins rather than searching the bamfiles for regions here because I want to ensure that each read is only mapped to a single feature

'''

trnasamplecounts  = defaultdict(int)
totalsamplecounts = defaultdict(int)

for currsample in samples:
    currbam = sampledata.getbam(currsample)
     
    try:
        #print >>sys.stderr, currbam
        pysam.index(""+currbam)
        bamfile = pysam.Samfile(""+currbam, "rb" )  
    except IOError as ( strerror):
        print >>sys.stderr, strerror
        sys.exit()
    
    for currread in getbamrange(bamfile, primaryonly = True):
        gotread = False
        totalsamplecounts[currsample] += 1
        for currbed in bedlist:
            for currfeat in featurelist[currbed].getbin(currread):
                if currfeat.coverage(currread) > 10:
                    counts[currsample][currbed] += 1
                    gotread = True
                    break
                    #print >>sys.stderr, currbam +":"+ currbed
        if gotread: 
            continue
        for currbed in locilist:
            for currfeat in trnaloci[currbed].getbin(currread):
                if currfeat.coverage(currread) > 10:
                    trnalocuscounts[currsample][currbed] += 1
                    #print >>sys.stderr, currfeat.bedstring()
                    gotread = True
                    break
        if gotread: 
            continue
        for currbed in trnalist:
            for currfeat in trnalist[currbed].getbin(currread):
                if currfeat.coverage(currread) > 10:
                    trnasamplecounts[currsample] += 1
                    trnacounts[currsample][currbed] += 1
                    fragtype = getfragtype(currfeat, currread)
                    if fragtype == "Whole":
                        
                        trnawholecounts[currsample][currbed] += 1
                    elif fragtype == "Fiveprime":
                        trnafivecounts[currsample][currbed] += 1
                    elif fragtype == "Threeprime":
                        trnathreecounts[currsample][currbed] += 1
                    elif fragtype == "Trailer":
                        trnatrailercounts[currsample][currbed] += 1
                    gotread = True
                    break
                        #print >>sys.stderr, str(currread.start - currfeat.start)+"-"+str(currread.end - currfeat.start)  
                        #print >>sys.stderr, str(currfeat.start - currfeat.start)+"-"+str(currfeat.end - currfeat.start)
                        #print >>sys.stderr, "****"
        if gotread: 
            continue
        othercounts[currsample] += 1
    
def sumsamples(countdict,sampledata, repname, currfeat = None, sizefactors = defaultdict(lambda: 1)):
    if currfeat is None: #To account for the "other" counts, which don't have a feature
        return sum(countdict[currsample]/sizefactors[currsample] for currsample in sampledata.getrepsamples(repname))
    else:
        return sum(countdict[currsample][currfeat]/sizefactors[currsample] for currsample in sampledata.getrepsamples(repname))
#tRNA-Ser-AGA-1-1    
    
if combinereps:
    replicates = list(sampledata.allreplicates())
    print "\t".join(replicates)
    for currbed in trnalist:
        
        if countfrags:
            #sumsamples(trnafivecounts,sampledata,currrep, currfeat)
            
            print "tRNA_wholecounts\t"+"\t".join(str(sumsamples(trnawholecounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
            print "tRNA_fiveprime\t"+"\t".join(str(sumsamples(trnafivecounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
            print "tRNA_threeprime\t"+"\t".join(str(sumsamples(trnathreecounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
            print "tRNA_other\t"+"\t".join(str(sumsamples(trnacounts,sampledata,currrep, currbed, sizefactors = sizefactor) - (sumsamples(trnafivecounts,sampledata,currrep, currbed, sizefactors = sizefactor) + sumsamples(trnathreecounts,sampledata,currrep, currbed, sizefactors = sizefactor) + sumsamples(trnawholecounts,sampledata,currrep, currbed, sizefactors = sizefactor))) for currrep in replicates)
        else:
            
            print currbed+"\t"+"\t".join(str(sumsamples(trnacounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)        
        
    
    for currbed in locilist:
        #print >>sys.stderr, currbed 
        print currbed+"\t"+"\t".join(str(sumsamples(trnalocuscounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
        
    for currbed in bedlist:  
        print currbed+"\t"+"\t".join(str(sumsamples(counts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
    print "other"+"\t"+"\t".join(str(sumsamples(othercounts,sampledata,currrep, sizefactors = sizefactor)) for currrep in replicates)
else:
    print "\t".join(samples)
    for currbed in trnalist:
        
        if countfrags:
            print "tRNA_wholecounts\t"+"\t".join(str(trnawholecounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
            print "tRNA_fiveprime\t"+"\t".join(str(trnafivecounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
            print "tRNA_threeprime\t"+"\t".join(str(trnathreecounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
            print "tRNA_other\t"+"\t".join(str((trnacounts[currsample][currbed] - (trnathreecounts[currsample][currbed] + trnafivecounts[currsample][currbed] + trnawholecounts[currsample][currbed]))/sizefactor[currsample]) for currsample in samples)
        else:
            print currbed+"\t"+"\t".join(str(trnacounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
        
        
    
    for currbed in locilist:
        print os.path.basename(currbed)+"\t"+"\t".join(str(trnalocuscounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
        
    for currbed in bedlist:
        print os.path.basename(currbed)+"\t"+"\t".join(str(counts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
    print "other"+"\t"+"\t".join(str(othercounts[currsample]/sizefactor[currsample]) for currsample in samples)



    
if args.trnanormfile is not None:
    #samples trnasamplecounts.keys()
    trnanormfile = open(args.trnanormfile, "w")
    mean = 1.*sum(trnasamplecounts.values())/len(trnasamplecounts.values())
    print >>trnanormfile, "\t".join(samples)
    print >>trnanormfile, "\t".join(str(trnasamplecounts[currsample]/mean) for currsample in samples)
    
if args.allreadsnormfile is not None:    
    allreadsnormfile = open(args.allreadsnormfile, "w")
    mean = 1.*sum(totalsamplecounts.values())/len(totalsamplecounts.values())
    print >>allreadsnormfile,"\t".join(samples)
    print >>allreadsnormfile,"\t".join(str(totalsamplecounts[currsample]/mean) for currsample in samples)
