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
~/pythonsource/trnaseq/getcoverage.py --samplefile=agingshort.txt  --bedfile=sacCer3-maturetRNAs.bed --stkfile=sacCer3-align.stk 
'''

gapchars = set("-._~")
class readcoverage:
    def __init__(self, region):
        self.region = region
        self.samplereads = 0
        self.coverage = list()
        self.length = region.length()
        self.totalreads = 0
        for i in range(0,region.length()):
            self.coverage.append(0)
    def addread(self, read):
        self.totalreads += 1
        if self.region.strand == "+":
            start = max([0, read.start - self.region.start])
            end = min([self.length, read.end - self.region.start])
            
        else:
            start = max([0, self.region.end - read.end])
            end = min([self.length, self.region.end - read.start])
        
        for currpos in range(self.length):
            if start <= currpos <= end - 1:
                self.coverage[currpos] += 1
                
    def coveragelist(self):
        return self.coverage
    def coveragealign(self, alignment, gapoutput = "NA",sizefactor = 1):
        if len(self.coverage) != len(string.translate(alignment, None, str(gapchars))):
            print >>sys.stderr, "Alignment length does not match bed length"            
        i = 0
        for curr in alignment:
            #print >>sys.stderr, curr
            if curr in gapchars:
                yield gapoutput
            else:
                yield self.coverage[i]/sizefactor
                i += 1

count = 0

positions = list([0,1,2,3,4,5,6,7,8,9,0,11,12,13,14,15,16,17,'-',18,19,20,'-','-',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e','e',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86])
def gettnanums(trnaalign, margin = 0):
    trnanum = list()
    currcount = 0
    enum = 1
    gapnum = 1
    for i in range(margin):
        trnanum.append('head'+str(margin - i))
    for i, struct in enumerate(trnaalign.consensus):
        if struct in  set("+=*"):
            #special case to account for differences between loci/transcripts
            if currcount == 0 and struct == '=':
                currcount = 1
            if positions[currcount] == 'e':
                trnanum.append('e'+str(gapnum))
                enum += 1
                currcount += 1
            elif positions[currcount] == '-':
                trnanum.append('gap'+str(gapnum))
                gapnum += 1
                currcount += 1
            else:
                trnanum.append(str(positions[currcount]))
                currcount += 1
        else:
            trnanum.append('gap'+str(gapnum))
            gapnum += 1
    for i in range(margin):
        trnanum.append('tail'+str(i+1))
    return trnanum

parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--bedfile',  nargs='+', default=list(),
                   help='bed file with mature tRNA features')
parser.add_argument('--samplefile',
                   help='Sample file in format')
parser.add_argument('--stkfile',
                   help='Stockholm file')
parser.add_argument('--sizefactors',
                   help='Optional file including size factors that will be used for normalization')
parser.add_argument('--combinereps', action="store_true", default=False,
                   help='Sum samples that are replicates')
parser.add_argument('--edgemargin', type=int, default=0,
                   help='margin to add to feature coordinates')
'''
parser.add_argument('--trnapositions', action="store_true", default=False,
                   help='Use tRNA positions')
'''

args = parser.parse_args()
edgemargin = int(args.edgemargin)

sampledata = samplefile(args.samplefile)
samples = sampledata.getsamples()

alltrnas = list()
#print >>sys.stderr, " ".join(bamlist)
#gettnanums
trnastk = list(readrnastk(open(args.stkfile, "r")))[0]

positionnums = gettnanums(trnastk, margin = edgemargin)
trnastk = trnastk.addmargin(edgemargin)

try:
    basetrnas = list()
    for currfile in args.bedfile:
        basetrnas.extend(list(currbed for currbed in readbed(currfile)))       
except IOError as e:
    print >>sys.stderr, e
    sys.exit()

trnalist = list(curr.addmargin(edgemargin) for curr in basetrnas)
#featurelist = list(readbed(sys.argv[1]))
#featurelist = list(readbed(sys.argv[1]))

#./countcomplete.py hg19-nontrnas.bed hg19-tRNAs.bed hg19-complete-tRNApad.fa




sizefactor = dict()
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
'''

cmbuild --enone --hand trnamature-euk.cm MaturetRNAs.stk

transcript
.(((.(.(.((..,....,..<<.<<.___..____..._>>>>,.<...<<<<.__..__.___....>>>.>>,,..............,,,.<<<<<._______>>>.>...>...))....)))))::::

loci:
((..(..(...(.....(..(.....,....,.<<<<____.___...._>>>>.........................,.<...<.<.<<...___.__.._................................................................................................._>>>.>.....>,,<<<<<<<____.>>>>>>>,..,<<<.<<._______...>>>..>..>....))....).)...))):

loop = [\.\_]+
openstem = [\.\(\<]+?
closestem = [\.\)\>]+?
inter = [\.\,]
'''
#transcriptalign = re.compile(r"\.[\.\(]+?[]") 

            
            

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
        
    for i, currfeat in enumerate(basetrnas):
        allcoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
        for currread in getbamrange(bamfile, currfeat):
            #print >>sys.stderr,  basetrnas[i].bedstring()
            #print >>sys.stderr,  currfeat.bedstring()
            #print >>sys.stderr,  "****"
            if  basetrnas[i].coverage(currread) > 10:
                allcoverages[currsample][trnalist[i].name].addread(currread)
            else:
                pass
                #print >>sys.stderr, "***"
#print >>sys.stderr, allcoverages['dmSCDd12_1']["tRNA-Phe-GAA-2"].totalreads/sizefactor['dmSCDd12_1']
#print >>sys.stderr, allcoverages['dmSCDd1_1']["tRNA-Phe-GAA-2"].totalreads/sizefactor['dmSCDd1_1']


print "Feature"+"\t"+"Sample"+"\t"+"\t".join(positionnums)
for currfeat in trnalist:
    totalreads = sum(allcoverages[currsample][currfeat.name].totalreads for currsample in samples)
    if totalreads < 100 * len(samples):
        continue
    if args.combinereps:

        replicates = sampledata.allreplicates()
        for currrep in replicates:
            print currfeat.name+"\t"+currrep+"\t"+"\t".join(str(curr) for curr in sumsamples(allcoverages,sampledata,currrep,currfeat,sizefactor))
        
    else:
        for currsample in samples:
            print currfeat.name+"\t"+currsample+"\t"+"\t".join(str(curr) for curr in allcoverages[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name],sizefactor = sizefactor[currsample]))
        
